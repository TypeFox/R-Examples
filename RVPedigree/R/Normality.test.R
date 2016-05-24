##' Test for normality of the trait
##'
##' This function is used to test whether the phenotype is distributed
##'     normally, based on the Shapiro-Wilk test.
##' @title Test for normality of the trait/phenotype
##' @inheritParams RVPedigree
##' @inheritParams pvalue.VCC2
##' @param plot (logical) If set to TRUE a histogram will be plotted
##'     of the phenotype residuals after adjusting for covariates
##'     (default: FALSE).
##' @return A list with the following elements:
##'     \itemize{
##'     \item \code{SW.pvalue}: Shapiro-Wilk p-value that indicates
##'     whether the phenotype is distributed normally.
##'     \item \code{resid.0}: the residuals after regressing the
##'     phenotype onto X.
##'     }
##' @author Karim Oualkacha
##' @export
Normality.test <- function(y=NULL,
                           X=NULL,
                           pedigree=NULL,
                           plot=FALSE)
{

    father.id <- pedigree[,3]

    is.founder <- (father.id==0)
    indices.founders <- c(1:nrow(X))[is.founder]

    resid.0 <- lm(as.vector(y) ~ X)$residuals
    resid.1 <- resid.0[indices.founders]

    SW.pvalue <- shapiro.test(resid.1)$p.value

    if (plot =="TRUE")
    {
        hist(resid.0,
             freq = NULL,
             main=NULL,
             xlab="Phenotype residuals after adjusting for covariates")
    }

    return(list(SW.pvalue=SW.pvalue, resid.0=resid.0))
}
