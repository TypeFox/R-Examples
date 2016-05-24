invertinfo <- function(mat, silent=TRUE, stopOnError=FALSE){
#Stable inversion of the information matrix to obtain a variance-covariance matrix; makes use of the fact that mat is a symmetrical and positive definite matrix
##############################
  if(stopOnError){
    result <- list(vcov=chol2inv(chol(mat)), error_message=NULL)
  }else{
    result <- list(vcov=matrix(NA, ncol=ncol(mat)), error_message=NULL)
    vcov <- try(chol2inv(chol(mat)), silent=silent)
    if(class(vcov)=="try-error"){
      result$error_message <- paste("Error in invertinfo(mat) : \n", vcov[[1]], sep="")
      result$vcov <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
    }else{
      result$vcov <- vcov
    }  
  }
  dimnames(result$vcov) <- dimnames(mat)
  return(result)
}
