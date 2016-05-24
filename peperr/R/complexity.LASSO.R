complexity.LASSO <- function(response, x, full.data, ...){
   require(penalized)
   prof <- profL1(response=response, penalized=x, trace=FALSE, ...)
   lambda <- prof$lambda[which.max(prof$cvl)]
   lambda
}