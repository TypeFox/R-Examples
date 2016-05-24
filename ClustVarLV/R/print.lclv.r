#' Print the LCLV results
#' 
#' @param x an object of class \code{lclv}
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso LCLV
#' 
#' @export
#' 
print.lclv =  function (x, ...) 
{
  if (!inherits(x, "lclv")) 
    stop("non convenient object")
  p = x$param$p
  n = x$param$n
  EXTu=x$param$EXTu
  EXTr=x$param$EXTr
  method =  x$param$method
  
  cat("\n")
  cat(paste("number of variables: ", p), sep = " ")
  cat("\n")
  cat(paste("number of observations: ", n), sep = " ")
  cat("\n")
  cat(paste("consolidation for K in c(",x$param$nmax,":2)",sep = ""))
  cat("\n")
  cat("\n")
  cat("$tarbres: results of the clustering")
  cat("\n")
  cat("$partitionK or [[K]]: results of partition into K clusters")
  cat("\n")
  cat("     $clusters: groups membership (before and after consolidation)")
  cat("\n")
  cat("     $compt: latent components of the clusters (after consolidation) defined according to the Xr variables")
  cat("\n")
  cat("     $compc: latent components of the clusters (after consolidation) defined according to the Xu variables")
  cat("\n")
  cat("     $loading_v: loadings of the external Xr variables (after consolidation)")
  cat("\n")
  cat("     $loading_u: loadings of the external Xu variables (after consolidation)")
  cat("\n")
  
}