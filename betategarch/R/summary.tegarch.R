summary.tegarch <-
function(object, verbose=FALSE, ...)
{
  if(verbose){
    out <- object[-1]
  }else{
    if(is.null(object$hessian)){
      out <- object[-c(1,3:7)]
    }else{
      out <- object[-c(1,3:8)]
    }
  }
  return(out)
}
