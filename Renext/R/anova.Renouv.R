##*****************************************************************************
##      DO NOT ROXYGENISE THIS!!!
## 
##' Compute an analysis of deviance table for two nested Renouv
##' objects
##'
##' Of special interest is the case when the distribution of the
##' exeedances used in \code{object} is exponential while
##' \code{object1} uses a two-parameters alternative in the GPD
##' family. We know then that the convergence to the asymptotic
##' distribution is slow, and a numerical approximation of the
##' true distribution of the test statistic is used when possible,
##' i.e. when the objects do not use MAX or OTS data and the number
##' of exceedances is between 8 and 500.
##' 
##' @title Compute an analysis of deviance table for two nested Renouv
##' objects
##' 
##' @param object A Renouv model as fitted with \code{\link{Renouv}}.
##' 
##' @param object1 A Renouv model such that \code{object} is nested in
##' \code{object1}.
##' 
##' @param trace Level of verbosity. The value \code{0} prints nothing.
##' 
##' @param ... Not used yet.
##'
##' @return
##' An object of class "anova" inheriting from class "data.frame".
##'
##' @note
##' The deviance of the models should can not be interpreted: only the
##' difference of the deviance us useful.
##'
##'
##' @examples
##' ## test using historical data
##' fit1Exp <- Renouv(Garonne,  distname.y = "exponential", plot = FALSE)
##' fit1GPD <- Renouv(Garonne, distname.y = "GPD", plot = FALSE)
##' anova(fit1Exp, fit1GPD)
##' 
##' ## test without using historical data
##' x <- Garonne$OTdata$Flow
##' dur <- Garonne$OTinfo$effDuration
##'
##' fit2Exp <- Renouv(x,  threshold = 2700,  effDuration = dur,
##'                   distname.y = "exponential", plot = FALSE)
##' fit2GPD <- Renouv(x, threshold = 2700, effDuration = dur,
##'                   distname.y = "GPD", plot = FALSE)
##' anova(fit2Exp, fit2GPD)
##' 
anova.Renouv <- function(object, object1, trace = 1L, ...) {
    
    if (missing(object) || missing(object1))
        stop("two models must be specified in 'object' and 'object1'")
    
    model0 <- deparse(substitute(object))
    model1 <- deparse(substitute(object1))
    models <- list(model0, model1)
    
    if (trace) {
        cat("Models: \n")
        cat(sprintf("  o \'%s\' with exceedances dist. \"%s\"\n", model0,
                    object$distname.y))
        cat(sprintf("  o \'%s\' with exceedances dist. \"%s\"\n\n", model1,
                    object1$distname.y))
    }
    
    ## check that the models have the same data. Could be improved
    ## by adding more tests (data values, ...)
    for (elt in list("threshold", "nb.OT", "nobs",
                     c("history.MAX", "flag"), c("history.OTS", "flag"))) {
        if (object[[elt]] != object1[[elt]]) 
            stop("'object' and 'object1' must have the same '",
                 paste(elt, collapse = "$"), "' element") 
    }
    
    df <- c(object$df, object1$df)
    dev <- -2 * c(object$logLik, object1$logLik)
    dfDiff <- diff(df)
    devDiff <- diff(dev)
                  
    ## checks df
    if (dfDiff <= 0) stop("non-nested models. Bad order for models?")

    w <- -devDiff
    n <- object$nb.OT

    alternative <- tolower(object1$distname.y)

    ## Use 'LRExp.test' for exponential H0 and GPD alternative (no
    ## equality restriction should have been given on the GPD
    ## parameters)
    
    if ((object1$distname.y %in% c("gpd", "GPD", "lomax", "maxlo")) &&
        (object1$df == 3L) && (object$distname.y %in% c("exp", "exponential")))  {
        
        H0 <- "exponential"
        
        if (object$df < 2L) {
            stop("the exponential parameter can not be fixed",
                 " in the 'exp | gpd' case")
        }
        
        if (!object$history.MAX$flag && !object$history.OTS$flag) {
            
            if (n <= 500) {
                meth <- "num"
                method <- "numerical approximation"
            } else {
                meth <- "asymp"
                method <- "asymptotic approximation"
            }
            
        } else {
            ## the models embed MAX or OTS data.
            if (object$nobs < 50L) {
                warning("asymptotic distribution may not be accurate ",
                        "here")
            }
            meth <- "asymp"
            method <- "asymptotic approximation"   
        }
        
        Test <- LRExp.test(object$y.OT, alternative = alternative,
                           method = meth)
        pVal <- Test$p.value
          
    } else {
        method <- "asymptotic approximation"
        pVal <- pchisq(w, df = 1, lower.tail = FALSE) 
    }

    if (trace) cat("Method used: ", method, "\n\n")
    
    table <- data.frame(df, dev,  c(NA, w), c(NA, pVal))
    dimnames(table) <- list(models, c("df", "deviance",  "W", "Pr(>W)"))
    
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
    
    ## list(W = w, p.value = pVal, method = method)
    
}


