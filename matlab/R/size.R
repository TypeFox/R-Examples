###
### $Id: size.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Array dimensions.
###


##-----------------------------------------------------------------------------
setGeneric("size",
           function(X, dimen) {
               #cat("generic", match.call()[[1]], "\n")
               standardGeneric("size")
           })

setMethod("size",
          signature(X     = "vector",
                    dimen = "missing"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(vector, missing)", "\n")
#              return(as.size_t(length(X)))
# :NOTE: Incompatible with previous implementation but consistent with MATLAB
              callGeneric(matrix(X, nrow = 1))
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "missing"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(matrix, missing)", "\n")
              return(as.size_t(dim(X)))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "missing"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(array, missing)", "\n")
              return(as.size_t(dim(X)))
          })

setMethod("size",
          signature(X     = "vector",
                    dimen = "numeric"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(vector, numeric)", "\n")
              callGeneric(matrix(X, nrow = 1), dimen)
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "numeric"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(matrix, numeric)", "\n")
              callGeneric(X, as.integer(dimen))
          })

setMethod("size",
          signature(X     = "matrix",
                    dimen = "integer"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(matrix, integer)", "\n")
              return(getLengthOfDimension(X, dimen))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "numeric"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(array, numeric)", "\n")
              callGeneric(X, as.integer(dimen))
          })

setMethod("size",
          signature(X     = "array",
                    dimen = "integer"),
          function(X, dimen) {
              #cat(match.call()[[1]],
              #    "(", data.class(X), ", ", data.class(dimen), ")", "\n")
              return(getLengthOfDimension(X, dimen))
          })

setMethod("size",
          signature(X     = "missing",
                    dimen = "ANY"),
          function(X, dimen) {
              #cat(match.call()[[1]], "(missing, ANY)", "\n")
              stop(sprintf("argument %s missing", sQuote("X")))
          })

##-----------------------------------------------------------------------------
getLengthOfDimension <- function(X, dimen) {
    if (!is.array(X)) {
        stop(sprintf("argument %s must be matrix or array", sQuote("X")))
    }

    if (!(length(dimen) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("dimen")))
    } else if (dimen < 1) {
        stop(sprintf("argument %s must be a positive quantity",
                     sQuote("dimen")))
    }

    len <- if (dimen <= length(dim(X))) {
               dim(X)[dimen]
           } else {
               1	# singleton dimension
           }

    return(as.integer(len))
}

