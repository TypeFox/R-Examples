fit.LASSO <- function(response, x, cplx, ...){
   require(penalized)
   data <- as.data.frame(x)
   data$response <- response
   res <- penalized(response=response, data=data, lambda1=cplx, penalized=x, trace=FALSE, ...)
   attr(res, "response") <- response
   attr(res, "x") <- x
   res
}