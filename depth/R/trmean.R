trmean=function(x,alpha,W=function(dep,alpha){return(1)},method="Tukey",ndir=1000,approx=FALSE,eps=1e-8,...){
  
	if(is.data.frame(x)) x=as.matrix(x)
	if(is.list(x)) {
		m=length(x)
		n=length(x[[1]])
		y=matrix(0,n,m)
		for(i in 1:m){
			y[,i]=x[[i]]
			if(length(x[[i]])!=n){ stop("When using a list, each element must be a vector of the same length.") }
		}
		x=y
	}

 	match.arg(method,c("Tukey","Liu","Oja"))
 	W=match.fun(W)
	p=length(x[1,])
	n=length(x[,1])
	
	if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.\n")) }
 	if(p!=2 & method=="Liu"){ stop("Liu's depth can be calculated on bivaraite data sets only") }
	if(approx==TRUE & method!="Tukey"){ 
		warning("An approximate method is available only for Tukey's depth. Argument approx=TRUE ignored.") 
		approx=FALSE
	}
	
	dpth=apply(x,1,depth,x=x,method=method,ndir=ndir,approx=approx,eps=eps)
	
	keep=(dpth>=alpha)
	if(sum(keep)==0){
		warning("None of the data has the required depth.")
		return(NA)
	}
	
	y=x[keep,]
	wt=sapply(dpth,W,alpha=alpha,...)
	
	return(apply(rep(wt[keep],p)*y,2,sum)/sum(wt[keep]))

}

