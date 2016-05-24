#' Plot ANFIS training errors
#'
#' Plot the training error of the network. If trainingType is "on-line" then 
#' full pattern errors along the patterns of the whole training process; for 
#' a specific epoch or the epoch summary error.
#'
#' @param x ANFIS class object.
#' @param y not used but necessary for redefining the generic function.
#' @param epoch for on-line only: epoch == Inf the whole training error; 
#'  epoch == integer > 0 the give epoch trainings errors, epoch == 0 the abs 
#'  epoch training sum of errors.
#' @param ... plot additional parameters.
#'
#' @return output graphics.
#'
#' @include Anfis-plotMF.R
#' @exportMethod plot
#' @docType methods
#' @name plot
#' @rdname ANFIS-plot
#' @aliases plot,ANFIS-method
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="plot", signature="ANFIS", 
  definition=function(x, y, epoch=Inf, ...){
  ##Off-Line training
  invisible(switch(x@trainingType,
    "Not trained yet"=warning("Trying to plot a not trained ANFIS"),
    "trainHybridJangOffLine"=plot(getErrors(x), type="b", xlab="Epoch",
      ylab="Epoch Error", main="Training Hybrid Off-Line Jang Error", ...),
    "trainHybridOffLine"=plot(getErrors(x), type="b", xlab="Epoch",
      ylab="Epoch Error", main="Training Hybrid Off-Line Error", ...)
    )#switch
  )#invisible

  ##On-Line training
  if(x@trainingType == "trainHybridJangOnLine"){
    if(is.infinite(epoch)){
      ##The whole pattern training errors
      plot(x@errors, type="b", xlab="Pattern", ylab="Error", 
        main="Training Hybrid On-Line Jang Pattern errors", ...)
    }else{
      ##Training sum(abs(errors)) for each epoch
      if(epoch==0){
        epochError<-sapply(1:(length(x@errors)%/%nrow(x@X)), function(epoch){
          sum(abs(x@errors[1:nrow(x@X)+(epoch-1)*nrow(x@X)]))
        })
        plot(epochError, type="b", xlab="Epoch", ylab="Error", 
          main="Training Hybrid On-Line Jang Error", ...)
      }else{
        ##Training pattern errors for a given epoch
        if(epoch*nrow(x@X)<length(x@errors)){
          plot(x@errors[1:nrow(X)+nrow(X)*(epoch-1)], type="b", xlab="Epoch",
            ylab="Error", 
            main=paste("Training Hybrid Jang On-line error for patterns of",
            paste("epoch", epoch)), ...)
        }else{
          stop("Epoch out of range")
        }
      }
    }    
  }##On-Line training
})
