# @export
.ecdf <- function(x,complement=FALSE) {
  if(is.table(x)) {
    tabz <- x
  } else {
    tabz <- table(x)
  }
  if(complement) {
    yval <- cumsum(rev(as.numeric(tabz)))/sum(tabz)
    yval <- rev(yval)
  } else {
    yval <- cumsum(as.numeric(tabz))/sum(tabz)
  }
  xval <- as.numeric(names(tabz))
  
  return(cbind("x"=rev(xval),"cdf(x)"=rev(yval)))
}
