## Purpose: Entropy of a fitted "vlmc" object, see ?vlmc
## ------------------------------------------------------------------------
## $Id: entropy.R,v 1.12 2014/10/04 12:36:42 maechler Exp $
## Author: Martin Maechler, Date:  5 Apr 2000, 18:31

## Entropy  ===  Log[Likelihood] !
entropy <- function(object)
{
    if(!is.vlmc(object))
        stop("first argument must be a \"vlmc\" object; see ?vlmc")
    .Call(vlmc_entropy, object $ vlmc)
}

logLik.vlmc <- function(object, ...)
{
    r <- entropy(object)
    attr(r, "df") <- (object$alpha.len - 1) * unname(object$size[["context"]])
    attr(r, "nobs") <- object$n
    class(r) <- "logLik"
    r
}


### Maybe -- rather call this on 2 'vlmc' objects
entropy2 <- function(ivlmc1, ivlmc2, alpha.len = ivlmc1[1])
{
    ## Purpose: Entropy between two vlmc (sub) trees, see ?vlmc
    ## ------------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  10 Apr 2000

###-- untested -- maybe non-sense
    if(0 >= (alpha.len <- as.integer(alpha.len)))
        stop("alphabet length must be >= 1")
    if(ivlmc2[1] != alpha.len)
        stop("alphabet length differs for 2nd arg")

    .Call(vlmc_entropy2, ivlmc1, ivlmc2)
}
