int.simpson=function(fdataobj,equi=TRUE,method="TRAPZ"){
 if (!inherits(fdataobj, "fdata"))  fdataobj<-fdata(fdataobj)
 n<-nrow(fdataobj)
 out<-rep(NA,n)
 tt<-fdataobj$argvals
 if (equi & method=="TRAPZ"){   
  p<-length(tt)
#  w<-(c(0,tt[2:(p)-tt[1:(p-1)]])+c(tt[2:(p)-tt[1:(p-1)]],0))/(p-1)
  w<-(c(0,tt[2:(p)]-tt[1:(p-1)])+c(tt[2:(p)]-tt[1:(p-1)],0))/2
  out<-drop(fdataobj$data%*%w )
 }
 else{
 for (i in 1:n) {
   out[i]<-int.simpson2(tt,fdataobj$data[i,],equi=equi,method=method)
   }
   }
	return(out)
}


int.simpson2=function(x,y,equi=TRUE,method="TRAPZ"){
  n=length(x);ny=length(y)
	if (n!=ny) stop("Different length in the input data")
	if (n==2 || ny==2) method="TRAPZ"
  out <- switch(method,
        "TRAPZ" = {
   		idx=2:length(x)
	    value<-as.double((x[idx]-x[idx-1])%*%(y[idx]+y[idx-1]))/2
	  },"CSR" = {
     if (!equi){
        n=2*n-1
        app=approx(x,y,n=n);x=app$x;y=app$y}
     	h=(max(x)-min(x))/(n-1)
	   value=(h/3)*(y[n]+y[1]+2*sum(y[2*(1:((n-1)/2))+1])+4*sum(y[2*(1:((n-1)/2))]))
    },
    "ESR" = {
     if (!equi){
             n=2*n-1
             app=approx(x,y,n=n);x=app$x;y=app$y}
     h=(max(x)-min(x))/(n-1)
	   if (n<=4) stop("This method needs n>4")
     value=17*(y[1]+y[n])+59*(y[2]+y[n-1])+43*(y[3]+y[n-2])+49*(y[4]+y[n-3])
	   value=value+48*sum(y[5:(n-4)])
  	 value=(h/48)*value
    }
   )
	return(out)
}


