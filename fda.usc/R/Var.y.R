Var.y=function(y,S,Var.e=NULL){
  df=traza(S)
  if (is.null(Var.e)) {
     if (is.vector(y)) {
     	  n=length(y)
     	  y.est=S%*%y
        se=sum((y-y.est)^2,na.rm=TRUE)/(n-df)
	      Var.e=se*diag(ncol(S))
       }
     else {
     if (!is.fdata(y)) y<-fdata(y)
        y<-y[["data"]]
        n=ncol(S)
        y.est<-t(S%*%t(y))
	      Var.e<-t(y-y.est)%*%(y-y.est)/(n-df)
  	}
  }
  var.y=S%*%Var.e%*%t(S)
  return(var.y)
}
