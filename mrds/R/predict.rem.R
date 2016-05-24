#' @export
# see predict.ds for documentation
predict.rem <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,...){
  # Functions Used: predict.rem.fi and predict.ds 

  model <- object

  if(!is.null(newdata)){
    compute <- TRUE
    xmat <- newdata
  }else{
    xmat <- model$mr$data
  }

  xmat$distance <- 0
  ddfobj <- model$ds$ds$aux$ddfobj

  if(ddfobj$type=="gamma"){
    xmat$distance <- rep(apex.gamma(ddfobj),2)
  }

  xmat$offsetvalue <- 0
  p.0 <- predict(model$mr,newdata=xmat,integrate=FALSE,compute=compute)$fitted

  if(is.null(newdata)){
    pdot <- predict(model$ds,esw=FALSE,compute=compute,
                    int.range=int.range)$fitted
  }else{
    pdot <- predict(model$ds,newdata=newdata[newdata$observer==1,],
                    esw=FALSE,compute=compute,int.range=int.range)$fitted
  }

  fitted <- p.0*pdot

  if(is.null(newdata)){
    names(fitted) <- model$mr$mr$data$object[model$mr$mr$data$observer==1]
  }else{
    names(fitted) <- newdata$object[newdata$observer==1]
  }

  return(list(fitted=fitted))
}
