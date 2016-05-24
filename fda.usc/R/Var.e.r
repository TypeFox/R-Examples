Var.e=function(y,S){
  df=traza(S)
  if (is.vector(y)) {
	  n=length(y)
	  y.est=S%*%y
	  se=sum((y-y.est)^2,na.rm=TRUE)/(n-df)
	  var.e=se*diag(ncol(S))
  }
  else {
    if (!is.fdata(y)) y<-fdata(y)
    y<-y[["data"]]
	  n=ncol(S)
	  y.est<-t(S%*%t(y))
	  var.e<-t(y-y.est)%*%(y-y.est)/(n-df)
  	}
  return(var.e)
}
