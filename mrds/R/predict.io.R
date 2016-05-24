#' @export
# see predict.ds for documentation
predict.io <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,...){
  # Functions Used: predict.io.fi and predict.ds

  # change 6/7/05 jll: removed dopdot argument and took dopdot=FALSE
  # functionality and moved to plot.io which was the only place it was used.
  # Added newdata functionality in call to predict.ds.

  model <- object
  if(is.null(newdata)){
    xmat <- model$mr$mr$data
  }else{
    compute <- TRUE
    xmat <- newdata
  }

  xmat$distance <- 0
  ddfobj <- model$ds$ds$aux$ddfobj

  # for gamma models need to find where p(x)=1 (apex), set that as distance
  if(ddfobj$type=="gamma"){
    xmat$distance <- rep(apex.gamma(ddfobj),2)
  }

  # calculate ps for each part of the model
  xmat$offsetvalue <- 0
  p.0 <- predict(model$mr,newdata=xmat,integrate=FALSE,compute=compute)$fitted
  if(is.null(newdata)){
    pdot <- predict(model$ds,esw=FALSE,compute=compute,
                     int.range=int.range)$fitted
  }else{
    pdot <- predict(model$ds,newdata=newdata[newdata$observer==1,],
                    esw=FALSE,compute=compute,int.range=int.range)$fitted
  }

  # take the product for the overall p
  fitted <- p.0*pdot

  # set labels
  if(is.null(newdata)){
    names(fitted) <- model$mr$mr$data$object[model$mr$mr$data$observer==1]
  }else{
    names(fitted) <- newdata$object[newdata$observer==1]
  }

  return(list(fitted=fitted))
}
