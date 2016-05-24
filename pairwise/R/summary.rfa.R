#' @export summary.rfa
#' @title S3 Summary for rasch factor analysis
#' @description S3 summary method for object of class\code{"pair"}
#' @param object object of class\code{"pair"}
#' @param sortdif logical with default \code{sortdif=FALSE} wether to order items / persons by difficulty / ability.
#' @param ... other parameters passed trough
########################### hier die summary method fuer rfa #############################
summary.rfa<-function(object, sortdif=FALSE, ...){
  ### item analysis performed
  if(object$transposed==FALSE){
    cat("Principal Components over items of Rasch Residuals: \n") 
    if(sortdif==TRUE){
    threshold <- object$pers_obj$pair$threshold
    sigma <- object$pers_obj$pair$sigma
    #####
    threshold <- threshold[order(sigma), ]
    #####
    object$pca$loadings <- object$pca$loadings[order(sigma), ]
    
    sigma <- sort(sigma)
    cat("(ordered by item location) \n")
    }
  }
  ### person analysis performed
  if(object$transposed==TRUE){  
    cat("Principal Components over persons of Rasch Residuals: \n") 
    if(sortdif==TRUE){
      WLE <- object$pers_obj$pers$WLE
      #####
      object$pca$loadings <- object$pca$loadings[order(WLE), ]
      
      WLE <- sort(WLE)
      cat("(ordered by person ability) \n")
    }
  }
  return(object$pca)
  #print(round(data.frame(object[2:1]),5))
}