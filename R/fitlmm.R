#' \code{Fit MLM to get the BLUP values}
#'
#' Get best linear unbiased predictors for phenotypic values.
#' Require package "lme4".
#'
#'
#' @param data Phenotypic data [data.frame], see example: "data/pheno.csv".
#' @param model [model], find help from lme4 function lmer.
#' @param which.factor The factor used for BLUP calculation, normally Genotype or Line [character].
#' @param outfile Output file name [character="blup.csv"].
#'
#' @return return fitted model.
#'
#' @examples
#' pheno <-  read.csv("data/pheno.csv", header=T)
#' fit <- get_BLUP(data = df, model = Brix ~ (1|Line) + (1|Loc) + (1|Year) + (1|Line:Loc) + (1|Line:Year),
#'                  which.factor = "Line", outfile="data/blup.csv")
#'
#' @export
get_BLUP <- function(data = df, model = Brix ~ (1|Line) + (1|Loc) + (1|Year) + (1|Line:Loc) + (1|Line:Year),
                     which.factor = "Line", outfile="data/blup.csv") {

    fit <- lmer(model, data=data)

    myblup <-  ranef(fit)
    # look at output structure
    #str(brixblup)
    # extract blup for line
    genoblup <- as.data.frame(myblup[which.factor])
    names(genoblup)[1] <- "value"

    feff <- fixef(fit)
    genoblup$value <- genoblup$value+feff
    # see the structure of the blup for each line
    #str(brixlineblup)
    # save the brixlineblup output to a separate .csv file
    write.csv(genoblup,  file=outfile)

    return(fit)
}

#' \code{Get Broad sense heritability from BLUP model.}
#'
#' Run the model fit using function "get_BLUP" and get H2.
#'
#'
#' @param fit Model fit returned from get_BLUP [model].
#' @param numerator Numerator of Vg/Ve, [char].
#' @param denominator Denominator of Vg/Ve [data.table], 2nd, degree of freedom.
#'        i.e. data.frame(f=c("Line", "Line:Year", "Line:Loc", "Residual"), df=c(1, 2, 2, 4))
#'
#' @return return A value of heritability.
#'
#' @examples
#' pheno <-  read.csv("data/pheno.csv", header=T)
#' fit <- get_BLUP(data = df, model = Brix ~ (1|Line) + (1|Loc) + (1|Year) + (1|Line:Loc) + (1|Line:Year),
#'                  which.factor = "Line", outfile="data/blup.csv")
#' get_H2(fit, numerator="Line",
#'        denominator=data.frame(f=c("Line", "Line:Year", "Line:Loc", "Residual"),
#'        df=c(1, 2, 2, 4)))
#'
#' @export
get_H2 <- function(fit, numerator="Line",
                   denominator=data.frame(f=c("Line", "Line:Year", "Line:Loc", "Residual"),
                                          df=c(1, 2, 2, 4))){
    df <- as.data.frame(VarCorr(fit))
    vg <- df[df$grp == numerator, "vcov"]

    ve <- 0
    for(i in 1:nrow(denominator)){
        v <- df[df$grp == denominator$f[i], "vcov"]/denominator$df[i]
        ve <- ve + v
    }
    return(vg/ve)
}
