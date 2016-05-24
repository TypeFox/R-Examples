## 2015-10-27
## not published; needs to copy many more arguments from 'object'
sighting.chat <- function (object, nsim = 1000) {
    if (!inherits(object, 'secr'))
        stop ("requires fitted secr model")
    CH <- object$capthist
    if (is.null(Tu(CH)) & is.null(Tm(CH)))
        stop("no unmarked or nonID sighting data in model")
    secr.fit (CH, mask = object$mask, fixed = object$fixed, details=list(nsim=nsim), 
              trace = FALSE, start = object)
}