##' Generates missing entries in data.frame/matrix
##'
##' @title Create random missing data
##' @param data data.frame
##' @param p Fraction of missing data in each column
##' @param cols Which columns (name or index) to alter
##' @param rowwise Should missing occur row-wise (either none or all selected columns are missing)
##' @param nafun (Optional) function to be applied on data.frame before return (e.g. \code{na.omit} to return complete-cases only)
##' @return data.frame
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
makemissing <- function(data,p=0.2,cols=seq_len(ncol(data)),rowwise=FALSE,nafun=function(x) x) {
  p <- rep(p,length.out=length(cols))
  if (!rowwise)
  for (i in seq_along(cols)) {
    data[rbinom(nrow(data),1,p[i])==1,cols[i]] <- NA
  }
  else
    data[which(rbinom(nrow(data),1,p)==1),cols] <- NA
  return(nafun(data))
}
