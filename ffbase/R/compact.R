#' Compact a ff vector or ffdf data frame
#'
#' Compact takes a ff vector and tries to use the smallest binary data type for this vector.
#' @aliases compact compact.ff compact.ffdf
#' @method compact ff
#' @method compact ffdf
#' @export
#' @param x \code{ff} or \code{ffdf} object
#' @param use.na \code{logical} if TRUE the resulting ff vector can contain NA, otherwise this is not checked
#' @param ... other parameters
#' @return compact cloned ff vector, or original if no compacting can be done
compact <- function(x, use.na=TRUE, ...){
   UseMethod("compact")
}

#' @export
#' @export compact.ff
compact.ff <- function(x, use.na=TRUE,...){
   vm <- which(.vmode == vmode(x))[1]
   if (vm > 9){
     return(x)
   }
   
   idx <- seq_len(vm)
   if (is.factor.ff(x)){
     r <- c(1, nlevels(x))
   } else {
     r <- range(x, na.rm=TRUE)
   }
   
   m <- (r[1] >= .vmin) & (r[2] <= .vmax)
   if (isTRUE(use.na)){
     m <- m & is.na(.vNA)
   }
   
   m <- which(m[idx])[1]
   if (m < vm){
     clone(x, vmode=.vmode[m])
   } else {
     x
   }
}

#' @export
#' @export compact.ffdf
compact.ffdf <- function(x, use.na=TRUE, ...){
   ret <- lapply(bit::physical(x), compact, use.na=use.na, ...)
   res <- do.call(ffdf, ret)
   close(x)
   res
}

# # testing 1,2,3
# irisf <- as.ffdf(iris)
# irisc <- compact(irisf)
# vmode(irisc)
