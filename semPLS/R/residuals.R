# 11.05.2011
residuals.sempls <- function(object, what=c("LVs", "MVs"), scale=c("scaled", "original"), total=FALSE){
    what <- match.arg(what)
    scale <- match.arg(scale)
    model <- object$model
    data <- object$data
    # LVs
    if(what=="LVs"){
        res <- object$factor_scores - predict(object, what=what, scale=scale, total=total)
    }
    # MVs
    else{
      if(scale=="scaled"){
        pdata <- predict(object, what=what, scale=scale, total=total)
        res <- data - pdata
      }
      else{
        data <- rescale(data)
        pdata <- predict(object, what=what, scale=scale, total=total)
        res <- data - pdata
      }
    }
    return(res)
}
