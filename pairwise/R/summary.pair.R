#' @export summary.pair
#' @title S3 Summary for item parameter
#' @description S3 summary method for object of class\code{"pair"}
#' @param object object of class\code{"pair"}
#' @param sortdif logical with default \code{sortdif=FALSE} wether to order items by difficulty.
#' @param ... other parameters passed trough
########################### hier die summary method fuer pair #############################
summary.pair<-function(object, sortdif=FALSE, ...){
  cat("Thurstonian thresholds and item location (sigma): \n") 
  if(sortdif==TRUE){
    threshold <- object$threshold
    #sb <- object$sb
    sigma <- object$sigma
    #####
    threshold <- threshold[order(sigma), ]
    #sb <- sb[order(sigma), ]
    sigma <- sort(sigma)
    object<-list(threshold=threshold,sigma=sigma)
    #class(object)<-c("PAIRl", "list")
    cat("(ordered by location) \n")
  }
  print(round(data.frame(object[2:1]),5))
}