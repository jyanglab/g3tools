#' \code{BayesCpi model for GS and GWAS}
#'
#' R codes based on Rohan Fernando and Dorian Garrick, Bayesian Methods Applied to GWAS.
#' All rights and credits go to Rohan Fernando and Dorian Garrick.
#' The code is not efficient for large scale genomic data.
#'
#'
#' @param genofile A data.frame contains snpid information, must contain col: mergeby. [data.frame, ["uid", ...]].
#' @param train_pheno A training data.frame contains columns of phenotypic traits and fixed effects.
#'                       Can be used for GWAS. [data.frame, ["uid", "pheno", ...]].
#' @param test_pheno A testing data.frame contains columns of phenotypic traits and fixed effects.
#'                      For GS. [data.frame, ["uid", "pheno", ...], =NULL].
#' @param mergeby The column id to merge the pheno and geno data. [chr, ="uid"]
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
#' @return return GWAS results.
#'
#' @examples
#'
#' @export
BayesCpi <- function(genofile = "data/bayes_geno.txt",
                    train_pheno = "data/bayes_train_pheno.txt",
                    test_pheno = NULL,
                    mergeby = "uid",
                    trait="TA",
                    seed = 12347,
                    chainLength = 1000,
                    probFixed = 0.999,
                    estimatePi = "yes",
                    dfEffectVar = 4,
                    nuRes = 4,
                    varGenotypic = 1,
                    varResidual = 1,
                    windowSize = 10,
                    outputFrequency = 100
                    ){

    # start from here
    library(data.table)
    set.seed(seed)
    pheno <- fread(train_pheno, header=TRUE)
    gplist <- read_geno_pheno(genofile, train_pheno, test_pheno, mergeby)

    gp <- gplist[[1]]
    UID = unname(as.matrix(gp[,1]))                   # First field is unique identifier
    y <- gp[, trait]                                  # trait values, trait of interest
    Z <- gp[, (ncol(pheno)+1): ncol(gp)]              # genotypes 0, 1, 2
    Z <- unname(as.matrix(Z))
    #Z        = unname(as.matrix((Z + 10)/10));       # Recode genotypes to 0, 1, 2 (number of B alleles)
    markerID = colnames(gp)[(ncol(pheno)+1):ncol(gp)]               # Remember the marker locus identifiers
    remove(gp)


    if(length(gplist) > 1){
        commonTestData <- gplist[[2]]
        testID = unname(as.matrix(commonTestData[,1]))                                 # First field is unique identifier
        yTest        = commonTestData[, trait]                                         # Second field is trait values
        ZTest        = commonTestData[, (ncol(pheno)+1) : ncol(commonTestData)]        # Remaining fields are GenSel-coded genotypes
        ZTest        = unname(as.matrix((ZTest)));                                     # Recode genotypes to 0, 1, 2 (number of B alleles)
        remove(commonTestData)
    }


    nmarkers = ncol(Z)                                              # number of markers
    nrecords = nrow(Z)                                              # number of animals

    # center the genotype matrix to accelerate mixing
    markerMeans = colMeans(Z)                             # compute the mean for each marker
    Z = t(t(Z) - markerMeans)                             # deviate covariate from its mean
    p = markerMeans/2.0                                   # compute frequency B allele for each marker
    mean2pq = mean(2*p*(1-p))                             # compute mean genotype variance

    varEffects  = varGenotypic/(nmarkers*(1-probFixed)*mean2pq)         # variance of locus effects is computed from genetic variance
                                                          #(e.g. Fernando et al., Acta Agriculturae Scand Section A, 2007; 57: 192-195)
    scaleVar    = varEffects*(dfEffectVar-2)/dfEffectVar; # scale factor for locus effects
    scaleRes    = varResidual*(nuRes-2)/nuRes             # scale factor for residual variance
    logPi       = log(probFixed)                          # compute these once since probFixed does not change!
    logPiComp   = log(1-probFixed)

    numberWindows = nmarkers/windowSize                   # number of genomic windows
    numberSamples = chainLength/outputFrequency           # number of samples of genetic variances


    alpha           = array(0.0, nmarkers) # reserve a vector to store sampled locus effects
    meanAlpha       = array(0.0, nmarkers) # reserve a vector to accumulate the posterior mean of locus effects
    modelFreq       = array(0.0, nmarkers) # reserve a vector to store model frequency
    mu              = mean(y,na.rm=TRUE)   # starting value for the location parameter
    meanMu          = 0                    # reserve a scalar to accumulate the posterior mean
    geneticVar      = array(0,numberSamples) # reserve a vector to store sampled genetic variances

    # reserve a matrix to store sampled proportion proportion of variance due to window
    windowVarProp   = matrix(0,nrow=numberSamples,ncol=numberWindows)
    sampleCount     = 0                    # initialize counter for number of samples of genetic variances
    piMean          = 0                    # initialize scalar to accumulat the sum of Pi
    delta           = array(1, nmarkers)   # vector to indicate locus effect in model or not

    # adjust y for the fixed effect (ie location parameter)
    y[is.na(y)] <- mu                      # impute missing pheno to mean
    ycorr = y - mu

    #ZPZ=t(Z)%*%Z
    #zpz=diag(ZPZ)

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

        # sample delta and then the effect for each locus
        nLoci = 0                             # Counter for number of loci fitted this iteraction

        # sample effect for each locus
        for (locus in 1:nmarkers){
            ycorr = ycorr + Z[,locus]*alpha[locus]                    #phenotypes are adjusted for all but this locus
            rhs   = t(Z[,locus]) %*% ycorr                            #rhs of MME adjusted for all but this locus
            zpz   = t(Z[,locus]) %*% Z[,locus]                        #OLS component of MME for this locus
            v0    = zpz*vare                                          #Var(rhs|delta=0)
            v1    = zpz^2*varEffects + zpz*vare                       #Var(rhs|delta=1)

            logDelta0 = -0.5*(log(v0) + rhs^2/v0) + logPi             #This locus not fitted
            logDelta1 = -0.5*(log(v1) + rhs^2/v1) + logPiComp         #this locus fitted
            probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1))         #near 0 if locus poor, 1 if locus very good
            u = runif(1)                                              #sample from uniform distribution
            if(is.na(probDelta1)){
                alpha[locus] = 0                                      #Sample the locus effect from prior
                delta[locus] = 0
            }else{
                if(u < probDelta1){                                       #Accept the sample with Pr(delta=1|ELSE)
                    nLoci  = nLoci + 1                                    #Increment a counter for loci fitted this iteration
                    mmeLhs = zpz + vare/varEffects                        #Form the coefficient matrix of MME
                    invLhs = 1.0/mmeLhs                                   #Invert the coefficient matrix
                    mean   = invLhs*rhs                                   #Solve the MME for locus effect
                    alpha[locus] = rnorm(1, mean, sqrt(invLhs*vare))      #Sample the locus effect from data
                    ycorr = ycorr - Z[, locus]*alpha[locus]               #Adjust the data for this locus effect
                    meanAlpha[locus] = meanAlpha[locus] + alpha[locus]    #Accumulate the sum for posterior mean
                    modelFreq[locus] = modelFreq[locus] + 1               #Accumulate counter for acceptance of this locus
                    delta[locus] = 1                                      #record that this locus was fitted
                }else{
                    alpha[locus] = 0                                      #Sample the locus effect from prior
                    delta[locus] = 0                                      #record thtat this locus not fitted in the model
                }
            }

        }

        # sample the common locus effect variance
        varEffects = ( scaleVar*dfEffectVar + sum(alpha^2) )/rchisq(1, dfEffectVar+nLoci)
        if(estimatePi == "yes"){
            # sample Pi
            aa = nmarkers - nLoci + 1
            bb = nLoci + 1
            pi = rbeta(1, aa, bb)                                     #Sample pi from full-conditional
            piMean = piMean + pi                                      #Accumulative sum for posterior mean of pi
            logPi  = log(pi)
            logPiComp = log(1-pi)
        }
        if(iter %% outputFrequency == 0){
            message(sprintf("#> [bayesCpi]: iteration = %5d number of loci in model = %5d ...", iter, nLoci))
            aHatTest  = ZTest %*% meanAlpha/iter                      #compute genomic breeding values in test data
            corr      = cor(aHatTest, yTest)                          #correlation of yhat and y obs
            regr      = corr*sqrt(var(yTest)/var(aHatTest))           #regress yTest on aHatTest
            RSquared  = corr*corr
            message(sprintf("#> [bayesCpi]: Corr=%8.5f, Regr=%8.5f, R2=%8.5f", corr, regr, RSquared))

            sampleCount = sampleCount + 1
            geneticVar[sampleCount] = var(Z%*%alpha)                  #variance of genetic values
            wEnd = 0
            for(window in 1:numberWindows){
                wStart = wEnd + 1                                     #start of current window
                wEnd   = wEnd + windowSize                            #end of current window
                if(wEnd > nmarkers){
                    wEnd = nmarkers                                   #last window may be smaller than windowSize
                }
                windowVarProp[sampleCount, window] = var(Z[,wStart:wEnd]%*%alpha[wStart*wEnd])/geneticVar[sampleCount]
            }
        }
    }

    ### works fine above!

    meanMu = meanMu/chainLength
    meanAlpha = meanAlpha/chainLength
    modelFreq = modelFreq/chainLength
    piMean    = piMean/chainLength

    nTestUID  = nrow(ZTest)
    for(uid in 1:nTestUID){
        cat(sprintf("%15s\t%15.7f\n"), testID[uid], aHatTest[uid])
    }

    for(locus in 1:nmarkers){
        cat(sprintf("%15s %15.7f %10.7f \n", markerID[locus], meanAlpha[locus], modelFreq[locus]))
    }

    proc.time()-ptime

}

#' @rdname BayesCpi
read_geno_pheno <- function(genofile = "data/bayes_geno.txt",
                            train_pheno = "data/bayes_train_pheno.txt",
                            test_pheno=NULL,
                            mergeby){
    # start from here
    #library(data.table)
    #set.seed(seed)

    geno <- fread(genofile, header=TRUE, data.table=FALSE)
    tp   <- fread(train_pheno, header=TRUE, data.table=FALSE)
    # training data
    gp  <- merge(tp, geno, by=mergeby) # Merge genotype and phenotype by their uid
    message(sprintf("#> [read_geno_pheno] TRAINING DATA: loading [%s] individuals and [%s] SNPs", nrow(gp), ncol(gp)-ncol(tp)))

    if(!is.null(test_pheno)){
        testp <- fread(test_pheno, header=TRUE, data.table=FALSE)
        # test data
        td <- merge(testp,  geno, by=mergeby)
        message(sprintf("#> [read_geno_pheno] TEST DATA: loading [%s] individuals and [%s] SNPs", nrow(td), ncol(td)-ncol(td)))
        rm(testp)
    }

    # Free up space:
    rm(list=c("geno", "tp"))
    if(!is.null(test_pheno)){
        return(list(gp, td))
    }else{
        return(list(gp))
    }
}
