#' @export summary.pairSE
#' @title S3 Summary for item parameter with standard errors
#' @description S3 summary method for object of class\code{"pairSE"}
#' @param object object of class\code{"pairSE"}
#' @param sortdif logical with default \code{sortdif=FALSE} wether to order items by difficulty.
#' @param ... other parameters passed trough
########################### hier die summary method für pairSE #############################
summary.pairSE<-function(object, sortdif=FALSE, ...){
  cat("Thurstonian thresholds and item location (sigma) and their std. errors: \n")
  
#   if(sortdif==TRUE){
#     sorter <- order(object$parameter[,"sigma"])
#     object$parameter <- object$parameter[sorter,]
#     object$SE <- object$SE[sorter,]
#     cat("(ordered by location) \n")
#   }
  
  if(sortdif==TRUE){
    #sorter <- order(x$parameter[,"sigma"])
    sorter <- order(object$sigma)
    # x$parameter <- x$parameter[sorter,]
    object$SE <- object$SE[sorter,]
    object$SEsigma <- object$SEsigma[sorter]
    object$threshold <- object$threshold[sorter,] 
    object$sigma <- object$sigma[sorter] 
    cat("(ordered by location) \n")
  }
  
  
  print( t(data.frame(object))[  unlist(lapply((c(1:(dim(t(data.frame(object)))[1]/2))),function(x){seq(from=x,by=(dim(t(data.frame(object)))[1]/2),length.out=2)}))  ,]   )
}
