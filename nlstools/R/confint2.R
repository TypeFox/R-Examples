"confint2" <- function(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...)
{
    method <- match.arg(method)

    format.perc <- function(probs, digits)
      ## Not yet exported, maybe useful in other contexts:
      ## quantile.default() sometimes uses a version of it
      paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
            "%")
    
    ## Taken from confint.nls
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- seq_along(pnames)
    if (is.numeric(parm)) 
        parm <- pnames[parm]
    
    ## Taken from confint.default and modified slightly to use t-distribution
    asCI <- function(object, parm, level)
    {
        a <- (1 - level)/2
        a <- c(a, 1 - a)
#        pct <- stats:::format.perc(a, 3)
        pct <- format.perc(a, 3)
        fac <- qt(a, df.residual(object))
        
        parmInd <- match(parm, pnames)
        ci <- array(NA, dim = c(length(parmInd), 2), dimnames = list(parm, pct))
        ses <- sqrt(diag(vcov(object)))[parmInd]
        ci[] <- cf[parmInd] + ses %o% fac
        ci
    }

    ## Taken from confint.nls
    asProf <- function(object, parm, level)
    {
        message("Waiting for profiling to be done...")
        utils::flush.console()
        object <- profile(object, which = parm, alphamax = (1 - level)/4)
        confint(object, parm = parm, level = level, ...)    
    }

    switch(method, asymptotic = asCI(object, parm, level), profile = asProf(object, parm, level))
}

