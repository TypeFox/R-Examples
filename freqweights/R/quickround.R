## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-28 Mon 18:29 emilio on emilio-despacho>
## ============================================================

##' A wrapper of \code{\link[plyr]{round_any}} to round data sets to multiple of any number.
##'
##' The \code{x} may contain non numerical varibles (factor, character, logical). They will remain unchanged. The numerical or time variables will be rounded to a multiple fo the accuracy.
##' 
##' If \code{accuracy} is of length 1, then this value is applied to all the columns of the data set. Otherwise, its length must be the same as the number of columns of \code{x}, including non numerical variables. If any value if \code{NA}, the corresponding variable will remain unchanged.
##' @title Round data sets
##' @param x a \code{tbl} object or a numeric or POSIXct matrix
##' @param accuracy number to round to; for POSIXct objects, a number of seconds
##' @return A \code{tbl} object.
##' @seealso \code{\link[plyr]{round_any}}
##' @rdname quickround
##' @importFrom plyr llply round_any
##' @import dplyr
##' @export
##' @examples
##' quickround(iris,0.2)
##' quickround(iris[,1:3],c(0.2,0.5,1.0))
##'
##' tfq <- tablefreq(iris, vars=c("Sepal.Length","Species"))
##' a <- quickround(tfq, c(0.3, NA, NA))
##' b <- tablefreq(a, freq="freq")
##' b
quickround <- function(x, accuracy){
    if(any(class(x) %in% c('matrix', 'double', 'numeric')))
      x <- tbl_df(as.data.frame(x))
    nms <- tbl_vars(x)
    if(length(accuracy) == 1) {
      accuracy <- rep(accuracy,length(nms))
    }
       a <- llply(seq_len(length(nms)),
                     function(i){
                       if(!is.na(accuracy[i])) {
                         return(paste( nms[i] ," = myround(",nms[i],
                                      ", accuracy= ", accuracy[i],")"))
                       } else{
                           return(paste( nms[i] ," = ",nms[i]))
                       }
                     })
    args <- paste(a, collapse=", ")
    x %>% evaldp(mutate,args)
  }

## This function should go inside quickround. But there is a bug in dplyr
myround <- function(x, accuracy, f=round) {
      if(is.factor(x) || is.character(x) || is.logical(x)) {
        return(x)
      } else {
        return(round_any(x,accuracy,f))
      }
    }

