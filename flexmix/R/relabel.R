#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id$
#

setGeneric("dorelabel", function(object, perm, ...) standardGeneric("dorelabel"))

setMethod("dorelabel", signature(object="flexmix", perm="vector"), function(object, perm, ...) {
  object <- callNextMethod(object, perm)
  object@posterior$scaled <- object@posterior$scaled[,perm,drop=FALSE]
  object@posterior$unscaled <- object@posterior$unscaled[,perm,drop=FALSE]
  object@cluster <- order(perm)[object@cluster]
  object@size <- object@size[perm]
  names(object@size) <- seq_along(perm)
  object
})
        
setMethod("dorelabel", signature(object="FLXdist", perm="vector"), function(object, perm, ...) {
  if (length(perm) != object@k) stop("length of order argument does not match number of components")
  if (any(sort(perm) != seq_len(object@k))) stop("order argument not specified correctly")
  object@prior <- object@prior[perm]
  object@components <- object@components[perm]
  names(object@components) <- sapply(seq_along(object@components), function(k)
                              gsub("[0-9]+", k, names(object@components)[k]))
  object@concomitant <- dorelabel(object@concomitant, perm, ...)
  object
})

setMethod("dorelabel", signature(object="FLXP", perm="vector"), function(object, perm, ...) {
  object@coef <- object@coef[,perm,drop=FALSE]
  colnames(object@coef) <- sapply(seq_len(ncol(object@coef)), function(k)
                                  gsub("[0-9]+", k, colnames(object@coef)[k]))
  object
})

setMethod("dorelabel", signature(object="FLXPmultinom", perm="vector"), function(object, perm, ...) {
  object@coef <- object@coef[,perm,drop=FALSE]
  object@coef <- sweep(object@coef, 1, object@coef[,1], "-")
  colnames(object@coef) <- sapply(seq_len(ncol(object@coef)), function(k)
                                  gsub("[0-9]+", k, colnames(object@coef)[k]))
  object
})
          
setMethod("relabel", signature(object="FLXdist", by="character"),
function(object, by, which=NULL, ...)
{
    by <- match.arg(by, c("prior", "model", "concomitant"))

    if(by=="prior"){
        perm <- order(prior(object), ...)
    }
    else if(by %in% c("model", "concomitant")) {
        pars <- parameters(object, which = by)
        index <- grep(which, rownames(pars))
        if (length(index) != 1)
          stop("no suitable ordering variable given in 'which'")
        perm <- order(pars[index,], ...)
    }
    object <- dorelabel(object, perm=perm)
    object
})


setMethod("relabel", signature(object="FLXdist", by="missing"),
function(object, by, ...)
{   
    object <- relabel(object, by="prior", ...)
    object
})

setMethod("relabel", signature(object="FLXdist", by="integer"),
function(object, by, ...)
{
    if(!all(sort(by) == seq_len(object@k)))
        stop("if integer, ", sQuote("by"),
             " must be a permutation of the numbers 1 to ", object@k)

    object <- dorelabel(object, by)
    object
})
