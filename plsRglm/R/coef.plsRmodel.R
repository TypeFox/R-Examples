coef.plsRmodel <- function(object,type=c("scaled","original"),...)
{
  if(missing(type)){type="original"}
  if(type=="original"){
    res<-list(CoeffC=object$CoeffC,Coeffs=object$Coeffs)
  }
  if(type=="scaled"){
    res<-list(CoeffC=object$CoeffC,Std.Coeffs=object$Std.Coeffs)
  }
  class(res) <- "coef.plsRmodel"
  res
}
