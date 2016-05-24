#' Fast conditional mean
#' 
#' \code{byMean} works like a very fast version of tapply with (weighted) \code{FUN=mean} or \code{FUN=weighted.mean}.
#' @param x \code{numeric} vector to be averaged
#' @param by (list of) \code{factor(s)} for which the mean will be calculated
#' @param na.rm \code{logical} If \code{TRUE} \code{NA} values will be removed
#' @param weight \code{numeric} with of same length as \code{x}
#' @param ... not used
#' @return \code{array} with dimensions of \code{by}
#' @export
byMean <- function(x, by, na.rm=FALSE, weight=NULL, ...){
  
  if (!is.list(by)){
    index <- list(by)
    tdim <- sapply(index, nlevels)
    bin <- by
    nbins <- tdim
  } else {
  
    ndim <- length(index)
    tdim <- sapply(index, nlevels)
    
    ai <- do.call(cbind, index)
    dim(ai) <- c(nrow(ai), ndim)
    
    nbins <- prod(tdim)
    
    bin <- arrayIndex2vectorIndex(ai, tdim)
  }
  if (missing(weight) || is.null(weight)){
    bs <- .Call("binned_sum", as.numeric(x), as.integer(bin), as.integer(nbins), PACKAGE = "ffbase")
  } else {
    stopifnot(length(weight) == length(x))
    bs <- .Call("bySum", as.numeric(x), as.integer(bin), as.integer(nbins), as.numeric(weight), PACKAGE = "ffbase")    
  }
  
  res <- bs[,2] / bs[,1]
  is.na(res) <- bs[,1] == 0
  
  if (!na.rm){
    is.na(res) <- (bs[,3] > 0)
  }
  dim(res) <- c(tdim)
  dimnames(res) <- c(lapply(index, levels))
  res  
}


##### quick testing code ######
# x <- as.numeric(1:100000)
# bin <- as.factor(as.integer(runif(length(x), 1, 101)))
# x[1] <- NA
# weight <- runif(length(x))
# 
# byMean(1:10, factor(1:10), weight=10:1)
# 
# byMean(x, bin, na.rm=TRUE)
# 
# system.time({
#   replicate(50, {tapply(x, bin, function(i){c(mean=mean(i, na.rm=TRUE), na=sum(is.na(i)))})})
# })
# 
# 
# system.time({
#   replicate(50, {byMean(x, bin, na.rm=TRUE)})
# })
# 
# system.time({
#   replicate(50, {byMean(x, bin, na.rm=TRUE, weight=weight)})
# })
# 
# system.time({
#   replicate(50, {binned_sum(x, bin, nbins=100L)})
# })
# 
# require(aggregate)
# system.time({
#   replicate(50, {byMeans(x, bin, na.rm=TRUE, weight=weight)})
# })
# 
