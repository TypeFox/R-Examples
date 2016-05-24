#' @export
# see predict.ds for documentation
predict.trial <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,...){
  # Functions Used: predict.trial.fi and predict.ds
  model <- object

  if(is.null(newdata)){
    xmat <- model$data
    xmat <- xmat[xmat$observer==1&xmat$object %in%
                   as.numeric(names(model$fitted)),]
  }else{
    xmat <- newdata
    compute <- TRUE
  }

  if(!compute){
     return(list(fitted=model$fitted))
  }

  xmat$distance <- 0
  ddfobj <- model$ds$ds$aux$ddfobj

  if(ddfobj$type=="gamma"){
    xmat$distance <- as.vector(apex.gamma(ddfobj))
  }

  if(is.null(newdata)){
    pdot <- predict(model$ds,esw=FALSE,compute=compute,
                    int.range=int.range)$fitted
  }else{
    pdot <- predict(model$ds,newdata=newdata,esw=FALSE,compute=compute,
                    int.range=int.range)$fitted
  }

  fitted <- predict(model$mr$mr,newdata=xmat,type="response")*pdot

  if(is.null(newdata)){
    names(fitted) <- names(model$fitted)
  }else{
    names(fitted) <- newdata$object[newdata$observer==1]
  }

  return(list(fitted=fitted))
}
