detect.normality.outliers=function(x,alpha=0.05){
  x=x[!is.na(x)]
  mu.x=median(x,na.rm=TRUE)
  var.x=mad(x,na.rm=TRUE)
  x.ok=(1-pnorm(abs(x-mu.x),mean=0,sd=sqrt(var.x))^(length(x)))>alpha
  return(x.ok)
}
