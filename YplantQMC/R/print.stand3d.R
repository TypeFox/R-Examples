#'@method print stand3d
#'@S3method print stand3d
print.stand3d <- function(x,...){
  
  n <- length(x$plants)
  cat("Yplant - stand object (class \'stand3d\').\n\n")
  cat(paste(c(rep("-",30),"\n"),collapse=""))
  cat("Contains",n,"plants.\n")
  cat("Leaf area index :",round(x$LAI,3),"m2 m-2")
}
