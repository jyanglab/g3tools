#' \code{BayesCpi model for GS and GWAS}
#'
#' R codes based on Rohan Fernando and Dorian Garrick, Bayesian Methods Applied to GWAS.
#' All rights and credits go to Rohan Fernando and Dorian Garrick.
#' The code is not efficient for large scale genomic data.
#'
#'
#' @param genofile A data.frame contains snpid information, must contain col: uid. [data.frame, ["uid", ...]].
#' @param train_pheno A training data.frame contains columns of phenotypic traits and fixed effects.
#'                       Can be used for GWAS. [data.frame, ["uid", "pheno", ...]].
#' @param test_pheno A testing data.frame contains columns of phenotypic traits and fixed effects.
#'                      For GS. [data.frame, ["uid", "pheno", ...], =NULL].
#'
#' @param seed Seed for the random number generator. [integer, =12347].
#' @param chainLength Number of iterations for the MCMC chain. [interger, =11000].
#'
#' @param probFixed Parameter "pi" the probability SNP effect is zero [num, =.99].
#' @param estimatePi Wether to estimate pi or not. [chr, ="yes"].
#'
#' @param dfEffectVar Hyper parameter (degrees of freedom) for locus effect variance. [integer, =4].
#' @param nuRes Hyper parameter (degrees of freedom) for residual variance. [integer, =4].
#' @param varGenotypic Used to derive hyper parameter (scale) for locus effect variance. [num, =1].
#' @param varResidual Used to derive hyper parameter (scale) for residual variance. [num, =1].
#' @param windowSize Number of consecutive markers in a genomic window. [integer, =10].
#' @param outputFrequency Frequency for reporting performance and for computing genetic variances. [integer, =100]
#'
#'
#'
#' @return return A list of many values or a script to gcta command.
#'
#' @examples
#'
#' @export
simcrop <- function(genofile = "genotypes.dat",
                    train_pheno = "trainPhenotypes.dat",
                    test_pheno = NULL,

                    seed = 12347,
                    chainLength = 11000,
                    probFixed = 0.99,
                    estimatePi = "yes",
                    dfEffectVar = 4,
                    nuRes = 4,
                    varGenotypic = 1 ,
                    varResidual = 1,
                    windowSize = 10,
                    outputFrequency = 100
                    ){

    # start from here
    library(data.table)
    set.seed(seed)

    geno <- fread(genofile, header=TRUE, data.table=FALSE)
    tp   <- fread(train_pheno, header=TRUE, data.table=FALSE)
    gtp  <- merge(tp, geno, by="uid") # Merge genotype and phenotype by their uid

    if(!is.null(test_pheno)){
        testp    = fread(test_pheno, header=TRUE, data.table=FALSE)
        commonTestData       = merge(testp,  geno, by="uid")
        rm(testp)
    }

    # Free up space:
    rm(list=c("geno", "tp"))


    animalID = unname(as.matrix(commonTrainingData[,1]))                # First field is animal identifier
    y        = commonTrainingData[, 2]                                  # Second field is trait values
    Z        = commonTrainingData[, 3: ncol(commonTrainingData)]        # Remaining fields are GenSel-coded genotypes
    Z        = unname(as.matrix((Z + 10)/10));                          # Recode genotypes to 0, 1, 2 (number of B alleles)
    markerID = colnames(commonTrainingData)[3:ncol(commonTrainingData)] # Remember the marker locus identifiers
    remove(commonTrainingData)

    testID = unname(as.matrix(commonTestData[,1]))                  # First field is animal identifier
    yTest        = commonTestData[, 2]                              # Second field is trait values
    ZTest        = commonTestData[, 3: ncol(commonTestData)]        # Remaining fields are GenSel-coded genotypes
    ZTest        = unname(as.matrix((ZTest + 10)/10));              # Recode genotypes to 0, 1, 2 (number of B alleles)
    remove(commonTestData)

    nmarkers = ncol(Z)                                              # number of markers
    nrecords = nrow(Z)                                              # number of animals

    # center the genotype matrix to accelerate mixing
    markerMeans = colMeans(Z)                             # compute the mean for each marker
    Z = t(t(Z) - markerMeans)                             # deviate covariate from its mean
    p = markerMeans/2.0                                   # compute frequency B allele for each marker
    mean2pq = mean(2*p*(1-p))                             # compute mean genotype variance

    varEffects  = varGenotypic/(nmarkers*mean2pq)         # variance of locus effects is computed from genetic variance
    #(e.g. Fernando et al., Acta Agriculturae Scand Section A, 2007; 57: 192-195)
    scaleVar    = varEffects*(dfEffectVar-2)/dfEffectVar; # scale factor for locus effects
    scaleRes    = varResidual*(nuRes-2)/nuRes             # scale factor for residual variance

    numberWindows = nmarkers/windowSize                   # number of genomic windows
    numberSamples = chainLength/outputFrequency           # number of samples of genetic variances


    alpha           = array(0.0, nmarkers) # reserve a vector to store sampled locus effects
    meanAlpha       = array(0.0, nmarkers) # reserve a vector to accumulate the posterior mean of locus effects
    modelFreq       = array(0.0, nmarkers) # reserve a vector to store model frequency
    mu              = mean(y)              # starting value for the location parameter
    meanMu          = 0                    # reserve a scalar to accumulate the posterior mean
    geneticVar      = array(0,numberSamples) # reserve a vector to store sampled genetic variances
    # reserve a matrix to store sampled proportion proportion of variance due to window
    windowVarProp   = matrix(0,nrow=numberSamples,ncol=numberWindows)
    sampleCount     = 0                    # initialize counter for number of samples of genetic variances



    # adjust y for the fixed effect (ie location parameter)
    ycorr = y - mu


    ZPZ=t(Z)%*%Z
    zpz=diag(ZPZ)

    ptime=proc.time()
    # mcmc sampling
    for (iter in 1:chainLength){

        # sample residual variance
        vare = ( t(ycorr)%*%ycorr + nuRes*scaleRes )/rchisq(1,nrecords + nuRes)

        # sample intercept
        ycorr = ycorr + mu                    # Unadjust y for the previous sample of mu
        rhs    = sum(ycorr)                   # Form X'y
        invLhs = 1.0/nrecords                 # Form (X'X)-1
        mean = rhs*invLhs                     # Solve (X'X) mu = X'y
        mu = rnorm(1,mean,sqrt(invLhs*vare))  # Sample new location parameter
        ycorr = ycorr - mu                    # Adjust y for the new sample of mu
        meanMu = meanMu + mu                  # Accumulate the sum to compute posterior mean

        # sample effect for each locus
        for (locus in 1:nmarkers){

            rhs=t(Z[,locus])%*%ycorr +zpz[locus]*alpha[locus]
            mmeLhs = zpz[locus] + vare/varEffects                        # Form the coefficient matrix of MME
            invLhs = 1.0/mmeLhs                                   # Invert the coefficient matrix
            mean = invLhs*rhs                                     # Solve the MME for locus effect
            oldAlpha=alpha[locus]
            alpha[locus]= rnorm(1,mean,sqrt(invLhs*vare))         # Sample the locus effect from data
            ycorr = ycorr + Z[,locus]*(oldAlpha-alpha[locus]);               # Adjust the data for this locus effect
            meanAlpha[locus] = meanAlpha[locus] + alpha[locus];   # Accumulate the sum for posterior mean
        }

        # sample the common locus effect variance
        varEffects = ( scaleVar*dfEffectVar + sum(alpha^2) )/rchisq(1,dfEffectVar+nmarkers)

    }

    proc.time()-ptime

}
