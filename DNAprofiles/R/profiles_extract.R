#' @name profiles extract
#' @title profiles extract
#' @method [ profiles
#' @usage \method{[}{profiles}(x,i,j,...)
#' @param x profiles object
#' @param i vector or matrix
#' @param j (optionally) vector
#' @param ... passed on to \code{[.matrix}
#' @description Profiles are stored in a profiles object, which is merely an integer matrix together with allelic frequencies stored as an attribute "freqs".
#' @examples data(freqsNLsgmplus)
#' x <- sample.profiles(N=10,freqsNLsgmplus)
#' y <- unclass(x);attr(y,"freqs") <- NULL
#' 
#' stopifnot(identical(as.vector(x[1,]),as.vector(y[1,])))
#' stopifnot(identical(as.vector(x[1:10]),as.vector(y[1:10])))
#' stopifnot(identical(as.vector(x[5,5]),as.vector(y[5,5])))
#' stopifnot(identical(as.vector(x[,1:10]),as.vector(y[,1:10])))
#' stopifnot(identical(as.vector(x[,]),as.vector(y[,])))
#' stopifnot(identical(as.vector(x[cbind(2,1:10)]),as.vector(y[cbind(2,1:10)])))
#' @export
"[.profiles" <- function(x, i,j,...){  
  y <- NextMethod(.Generic)
  attr(y,"freqs") <- attr(x,"freqs")
  class(y) <- .Class
  y
}

# 
# "[.profiles" <- function(x, i,j,drop=TRUE){  
#   Narg <- nargs() - !missing(drop)  # number of arg from x,i,j
#   has.j <- !missing(j)
#   
#   if (Narg<3L){ # x[], x[1], x[1:2] or x[cbind(1,1:2)]
#     ret <- unclass(x)[i,drop=drop] # extract like a matrix    
#   }else{ # x[,], x[,j], etc.
#     ret <- unclass(x)[i,j,drop=drop] # extract like a matrix        
#   }
#   
#   class(ret) <- c("profiles",class(ret)) # return as a profiles
#   attr(ret,"freqs") <- attr(x,"freqs")
#   ret
# }