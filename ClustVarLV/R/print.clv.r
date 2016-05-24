#' Print the CLV results
#' 
#' @param x an object of class \code{clv}
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso CLV
#' 
#' @export
#' 
print.clv =  function (x, ...) 
{
  if (!inherits(x, "clv")) 
    stop("non convenient object")
  p = x$param$p
  n = x$param$n
  EXTu=x$param$EXTu
  EXTr=x$param$EXTr
  method =  x$param$method
  strategy = x$param$strategy
  if (is.null(x$param$nmax)) {
    cas="clvkmeans"
  } else {
    cas="clv"
  }
  
  if (cas=="clv") {
    cat("\n")
    cat(paste("number of variables: ", p), sep = " ")
    cat("\n")
    cat(paste("number of observations: ", n), sep = " ")
    cat("\n")
    if (method==1) cat("measure of proximity: squared covariance")
    if (method==2) cat("measure of proximity: covariance")
    cat("\n")
    cat(paste("consolidation for K in c(",x$param$nmax,":2)",sep = ""))
    cat("\n")
    cat("\n")
    cat("$tarbres: results of the clustering")
    cat("\n")
    cat("$partitionK or [[K]]: partition into K clusters")
    cat("\n")
    cat("    partition[[K]]$clusters: cluster's membership (1st line: before and 2nd line: after consolidation)")
    cat("\n")
    cat("    partition[[K]]$comp: latent components of the clusters (after consolidation),matrix of size (n x K)")
    
    if ((EXTu==1)|(EXTr==1)){
    cat("\n")
    cat("     $loading: loadings of the external variables (after consolidation)")
    }
    cat("\n")
  }
  
  if (cas=="clvkmeans") {
    cat("\n")
    cat(paste("number of variables: ", p), sep = " ")
    cat("\n")
    cat(paste("number of observations: ", n), sep = " ")
    cat("\n")
    if (method==1) cat("measure of proximity: squared covariance")
    if (method==2) cat("measure of proximity: covariance")
    cat("\n")
    cat(paste("number of clusters: ", x$param$K), sep = " ")
    cat("\n")
    if (strategy=="sparselv") cat(paste("number of var with loading null : ", length(which(unlist(x$sloading)==0))), sep = " ")
    if (strategy=="kplusone") cat(paste("number of var put aside : ", length(which(x$clusters[2,]==0))), sep = " ")
    cat("\n")
    cat("\n")
    cat("$tarbres: results of the clustering")
    cat("\n")
    cat("$clusters: cluster's membership (1st line: intial partition and 2nd line: partition at convergence")
    cat("\n")
    cat("$comp: latent components of the clusters;matrix of size (n x K)")
    if ((EXTu==1)|(EXTr==1)){
    cat("\n")
    cat("     $loading: loadings of the external variables (after consolidation)")
    }
    cat("\n")
  }
  
}