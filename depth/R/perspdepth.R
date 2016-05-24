perspdepth=function(x,method="Tukey",output=FALSE,tt=50,xlab="X",ylab="Y",zlab=NULL,col=NULL,...){

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

# Suppose n data points in p dimensions. 
# Plots depth of sample x.
# x= matrix n by p
# method= which depth to use


  require(rgl)
  match.arg(method,c("Tukey","Liu","Oja"))
  p=length(x[1,])
  n=length(x[,1])

  if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.")) }
  if(p!=2) { stop("Data must be bivariate.\n") }

  if(is.null(zlab)){ zlab=paste(method,"'s depth",sep="") }
 
  if(method=="Tukey"||method=="Liu"){
	y=x[,2]
	x=x[,1]
        minx=min(x)
        miny=min(y)
        ecx=max(x)-minx
	ecy=max(y)-miny
	xx=minx
	yy=miny

	for (i in 1:tt){
		xx=c(xx,minx+i/tt*ecx)
		yy=c(yy,miny+i/tt*ecy)
	}

	ans = .Fortran("iso3d",
			as.numeric(x),
			as.numeric(y),
			z=numeric((tt+1)^2),
			as.integer(n),
			as.integer(tt),
			as.integer(method=="Liu"),  # Controls depth defn:  0=Tukey , 1=Liu
			as.numeric(xx),
			as.numeric(yy),
			PACKAGE="depth")

	zz=matrix(ans$z,tt+1,tt+1)

      if(output==FALSE){
          if(is.null(col)){  col="lightblue" }
	  persp3d(xx,yy,zz,col=col,xlab=xlab,ylab=ylab,zlab=zlab,...) 
	  rep=NULL
	}
	if(output==TRUE){
		return(list(x=as.vector(xx),y=as.vector(yy),z=as.matrix(zz)))
	}
  
  }
   
  if(method=="Oja"){
        w=x
	y=x[,2]
	x=x[,1]
        minx=min(x)
        miny=min(y)
        ecx=max(x)-minx
	ecy=max(y)-miny
	xx=minx
	yy=miny

	for (i in 1:tt){
		xx=c(xx,minx+i/tt*ecx)
		yy=c(yy,miny+i/tt*ecy)
	}

	ans = .Fortran("ojaiso3d",
			as.numeric(w),
			z=numeric((tt+1)^2),
			as.integer(n),
			as.integer(tt),
			as.numeric(xx),
			as.numeric(yy),
			PACKAGE="depth")

	zz=matrix(ans$z,tt+1,tt+1)
 
        if(output==FALSE){ 
          if(is.null(col)){  col="lightblue" }
	  persp3d(xx,yy,zz,col=col,xlab=xlab,ylab=ylab,zlab=zlab,...) 
	  rep=NULL
	}
	if(output==TRUE){
		return(list(x=as.vector(xx),y=as.vector(yy),z=as.matrix(zz)))
	}
  }
  
  invisible()
}

