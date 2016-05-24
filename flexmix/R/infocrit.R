#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: infocrit.R 4922 2013-09-03 13:32:45Z gruen $
#

setMethod("nobs", signature(object="flexmix"),
function(object, ...) {          
  if (is.null(object@weights)) nrow(object@posterior$scaled) else  sum(object@weights)
})

setMethod("logLik", signature(object="flexmix"),
function(object, newdata, ...){
    if (missing(newdata)) {
        z <- object@logLik
        attr(z, "df") <- object@df
        attr(z, "nobs") <- nobs(object)
        class(z) <- "logLik"
    } else {
        z <- sum(log(rowSums(posterior(object, newdata = newdata, unscaled = TRUE))))
        attr(z, "df") <- object@df
        attr(z, "nobs") <- nrow(newdata)
        class(z) <- "logLik"
    }
    z
})

setMethod("ICL", signature(object="flexmix"),
function(object, ...){
    -2 * clogLik(object) + object@df * log(nobs(object))
})

setMethod("clogLik", signature(object="flexmix"),
function(object, ...){
    first <- if (length(object@group)) groupFirst(object@group) else TRUE
    post <- object@posterior$unscaled[first,,drop=FALSE]
    n <- nrow(post)
    sum(log(post[seq_len(n) + (clusters(object)[first] - 1)*n]))
})

setMethod("EIC", signature(object="flexmix"),
function(object, ...) {
    first <- if (length(object@group)) groupFirst(object@group) else TRUE
    post <- object@posterior$scaled[first,,drop=FALSE]
    n <- nrow(post)
    lpost <- log(post)
    if (any(is.infinite(lpost))) lpost[is.infinite(lpost)] <- -10^3
    1 + sum(post * lpost)/(n * log(object@k))
})
