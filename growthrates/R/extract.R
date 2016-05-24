#' Extract or Replace Parts of a 'multiple_fits' Object
#'
#' Operators to access parts of 'multiple_fits' objects
#'
#' @param x object of class multiple_fits
#' @param i numeric or character index
#' @param j NULL (for compatibility with other uses of  \code{[} or \code{[[})
#' @param drop	If \code{TRUE} the result is coerced to the lowest possible
#'   dimension
#' @param \dots optional arguments passed to \code{[}
#'
#' @examples
#'
#' data(bactgrowth)
#' L <- all_splines(value ~ time | strain + conc +replicate, data=bactgrowth)
#'
#' coef(L[[1]])
#'
#' plot(L[["R:0:2"]])
#'
#' par(mfrow=c(2, 2))
#' plot(L[1:4])
#'
#' @rdname extract
#' @exportMethod "["
#'
setMethod("[", c(x="multiple_fits", i="ANY", j="missing"),
          function(x, i, j=NULL, ...) {
            ## if (!is.null(j))
            ##   stop("incorrect number of subscripts")
            #if (length(i) == 1) {
            #  x@fits[i=i]
            #} else {
              new("multiple_fits",
                  fits = x@fits[i=i],
                  grouping = x@grouping
                  )
            #}
          })

#' @rdname extract
#' @exportMethod "[["
#'
setMethod("[[", c(x="multiple_fits", i="ANY", j="missing"),
         function(x, i, j, ...) {
           ## if (!is.null(j))
           ##   stop("incorrect number of subscripts")
           if (length(i) > 1)
             stop("[[ can only be used to select one single element")
            x@fits[[i]]
          })
