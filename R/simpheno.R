#' \code{Simulate QTL phenotype with normal distribution}
#'
#' Simulate phenotype using genotype data given number of QTLs and heritability.
#'
#'
#' @param geno genotype data, col1=genoid, col2 and after=snpid, coding: 0,1,2 (no missing data allowed) [data.frame].
#' @param h2 Broad sense heritability of the trait [numeric(0-1)].
#' @param alpha Alpha.
#' @param nqtl number of QTL [interger].
#' @param distribution [character=norm]
#'
#' @return return A list of many values.
#'
#' @examples
#' geno <- read.table("data/geno.txt", header=TRUE)
#' pheno <- sim_qtl_pheno(geno, h2=0.7, alpha=0.5, nqtl=10, distribution="norm")
#' y <- pheno[['y']]
#'
#' @export
sim_qtl_pheno <- function(geno, h2, alpha, nqtl, distribution="norm"){

    X <- geno[, -1]
    n <- nrow(X)
    m <- ncol(X)

    #Sampling QTN
    QTN.position <- sample(m, nqtl, replace=F)
    SNPQ <- as.matrix(X[, QTN.position])

    message(sprintf("[g3tools:simpheno], read in [ %s ] SNPs for [ %s ] plants and simulated [ %s ] QTLs",
                    m, n, nqtl))

    #Simulate phenotype
    if(distribution == "norm"){
        addeffect <- rnorm(nqtl,0,1)
    }else{
        addeffect <- alpha^(1:nqtl)
    }

    effect <- SNPQ%*%addeffect
    effectvar <- var(effect)
    residualvar <- (effectvar - h2*effectvar)/h2
    residual <- rnorm(n, 0, residualvar)
    y <- as.data.frame(effect + residual)

    return(list(addeffect = addeffect, y=y, add = effect, residual = residual,
                QTN.position=QTN.position,SNPQ=SNPQ))
}
