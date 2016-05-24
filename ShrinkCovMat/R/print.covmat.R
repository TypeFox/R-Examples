print.covmat <-
function (x, ...) 
{
 cat("SHRINKAGE ESTIMATION OF THE COVARIANCE MATRIX", "\n")
 cat("\nEstimated Optimal Shrinkage Intensity =",round(x$lambdahat,4),"\n")
 cat("\nEstimated Covariance Matrix [1:5,1:5] =\n")
 print(round(x$Sigmahat[1:min(5,nrow(x$Sigmahat)),1:min(5,nrow(x$Sigmahat))],4))
 cat("\nTarget Matrix [1:5,1:5] =\n")
 print(round(x$Target[1:min(5,nrow(x$Target)),1:min(5,nrow(x$Target))],4))
}
