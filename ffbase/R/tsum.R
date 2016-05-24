#' Fast conditional sum
#' 
#' \code{condSum} works like a very fast version of tapply with \code{FUN=sum}.
#' @param x \code{numeric} vector to be summed
#' @param index (list of) \code{factor(s)} for which the sum will be calculated
#' @param na.rm \code{logical} If \code{TRUE} \code{NA} values will be removed
#' @param ... not used
#' @return \code{array} with dimensions of \code{index}
#' @export
condSum <- function(x, index, na.rm=FALSE, ...){
   
   if (!is.list(index)){
     index <- list(index)
   }
   
   ndim <- length(index)
   tdim <- sapply(index, nlevels)
   nbins <- prod(tdim)
   
   ai <- do.call(cbind, index)
   dim(ai) <- c(nrow(ai), ndim)
   bin <- arrayIndex2vectorIndex(ai, tdim)
   
   bs <- binned_sum(x, bin, nbins)
   
   res <- bs[,2]
   if (!na.rm){
     is.na(res) <- (bs[,3] > 0)
   }
   
   # TODO add NA's
   dim(res) <- c(tdim)
   dimnames(res) <- c(lapply(index, levels))
   res
}

#' Fast conditional mean
#' 
#' \code{condMean} works like a very fast version of tapply with \code{FUN=mean}.
#' @param x \code{numeric} vector to be averaged
#' @param index (list of) \code{factor(s)} for which the mean will be calculated
#' @param na.rm \code{logical} If \code{TRUE} \code{NA} values will be removed
#' @param ... not used
#' @return \code{array} with dimensions of \code{index}
#' @export
condMean <- function(x, index, na.rm=FALSE, ...){
  
  if (!is.list(index)){
    index <- list(index)
  }
  
  ndim <- length(index)
  tdim <- sapply(index, nlevels)
  
  ai <- do.call(cbind, index)
  dim(ai) <- c(nrow(ai), ndim)
  
  nbins <- prod(tdim)
  
  bin <- arrayIndex2vectorIndex(ai, tdim)
  
  bs <- binned_sum(x, bin, nbins)
  
  res <- bs[,2] / bs[,1]
  is.na(res) <- bs[,1] == 0
  
  if (!na.rm){
    is.na(res) <- (bs[,3] > 0)
  }
  
  dim(res) <- c(tdim)
  dimnames(res) <- c(lapply(index, levels))
  res
}

# x <- c(1:11,NA)
# f <- factor(sample(c("a", "b"), 12, replace=TRUE))
# f2 <- factor(sample(c("c", "d"), 12, replace=TRUE))
# 
# condSum(x, list(f=f,f2=f2))
# condSum(x, list(f=f,f2=f2), na.rm=TRUE)
# condMean(x, list(f=f,f2=f2), na.rm=TRUE)
