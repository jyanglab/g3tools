#' \code{Simulate QTL phenotype with normal distribution}
#'
#' Simulate crop phenotype using real genotype data for given number of QTLs, heritability.
#' The SNP effect can be drawn from normal and gamma distributions.
#'
#'
#' @param snpdf A data.frame contains snpid information, must contain col: snpid. [data.frame].
#' @param h2 Broad sense heritability of the trait. [numeric(0-1)].
#' @param nqtl number of QTL. [interger].
#' @param distribution Distribution of the effects. [character=norm/gamma]
#' @param usegcta For large set of SNP data, use external software GCTA to rapid simulate it. [TRUE/FALSE].
#' @param gctapar If usegcta==TRUE, here is the parameters passed to the gcta command. [vector]
#'                1. bfile, PLINK format
#'                2. path for causal.snpfile
#'                3. replication number
#'                4. path for output file.
#'
#'
#' @return return A list of many values or a script to gcta command.
#'
#' @examples
#' geno <- read.table("data/geno.txt", header=TRUE)
#' simcrop <- function(snpid, h2, nqtl, distribution="norm", usegcta=TRUE,
#'                     gctapar=c("test", "snplist.txt", 1, "outfile.txt"))
#' y <- pheno[['y']]
#'
#' gcta64  --bfile test  --simu-qt  --simu-causal-loci causal.snplist  --simu-hsq 0.5
#' --simu-rep 3  --keep test.indi.list --out test
#'
#' @export
simcrop <- function(snpdf, h2, nqtl, distribution="norm", usegcta=TRUE,
                    gctapar=c("test", "snplist.txt", 1, "outfile.txt")){

    m <- nrow(snpdf)

    #Sampling QTN
    idx <- sample(m, nqtl, replace=F)
    df <- snpdf[idx, ]

    #message(sprintf("[g3tools:simpheno], read in [ %s ] SNPs for [ %s ] plants and simulated [ %s ] QTLs",
    #                m, n, nqtl))

    #Simulate phenotype
    if(distribution == "norm"){
        addeffect <- rnorm(nqtl,0,1)
    }else if(distribution == "gamma"){
        addeffect <- rgamma(nqtl,1)*sample(c(-1,1), size=nqtl, replace=TRUE)
    }else{
        stop("### try norm or gamma ! ")
    }

    if(usegcta){
        df$effect <- addeffect
        write.table(df[, c("snpid", "effect")], gctapar[2], sep="\t",
                    row.names=FALSE, quote=FALSE, col.names=FALSE)
        cmd <- paste("gcta64 --bfile", gctapar[1],
                     "--simu-qt",
                     "--simu-causal-loci", gctapar[2],
                     "--simu-hsq", h2,
                     "--simu-rep", gctapar[3],
                     "--out", gctapar[4])
        return(cmd)

    ### not well implemented below
    }else{
        effect <- SNPQ%*%addeffect
        effectvar <- var(effect)
        residualvar <- (effectvar - h2*effectvar)/h2
        residual <- rnorm(n, 0, residualvar)
        y <- as.data.frame(effect + residual)

        return(list(addeffect = addeffect, y=y, add = effect, residual = residual,
                    QTN.position=idx, SNPQ=SNPQ))
    }

}
