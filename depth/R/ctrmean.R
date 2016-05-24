
ctrmean=function(x,alpha,eps=1e-8,mustdith=FALSE,maxdith=50,dithfactor=10,factor=.8){

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

	p=length(x[1,])
	n=length(x[,1])
	depth=ceiling(n*alpha)

	if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.\n")) }
 	if(p!=2){ stop("Depth contours can only be calculated on bivaraite data.") }
	if(length(depth)!=1 | round(depth)!=depth){ stop("The argument depth must be a single integer") }
	if(depth<1 | depth>ceiling(n/2)){ stop(message=paste("Depth must be an integer between 1 and ",ceiling(n/2))) }

	ndpth=1
	y=x[,2]
	x=x[,1]
	maxnum=floor(4*n*sqrt(n)+1)
	

	zz=.Fortran("halfmed",
			as.numeric(x),
			as.numeric(y),
			as.integer(n),
			integer(1),
			numeric(2),
			xcont=numeric(n*(n-1)/2),
			ycont=numeric(n*(n-1)/2),
			ncont=integer(ndpth),
			as.integer(depth),
			as.integer(ndpth),
			as.integer(1),
			as.integer(maxnum),
			err=integer(1),
			as.numeric(eps),
			as.numeric(dithfactor),
			as.integer(maxdith),
			as.integer(mustdith),
			missing=integer(ndpth),
			as.numeric(factor),
			PACKAGE="depth")
			

	if (zz$err==-5) { warning("Ventilation was used on the data.") }
	if (zz$err==-1) {
 		if(mustdith){ 
			warning("Data are not in general position. Dithering was used.") 
		} else {
			stop("Data are not in general position.")
		}	
	}
	if (zz$err==-2) { stop(message=paste(maxdith," ventilation steps and the data still fail to be in general position.")) }
	if (zz$err==-4) { warning("A non critical numerical error occurred during the calculation of the depth contour.") }
	if (zz$err>0) { stop(message=paste("CRITICAL numerical error in calculating contour")) }
	if(zz$ncont==0){ warning("The requested contour does not exist.") }
	n=zz$ncont
	x=zz$xcont[c(1:n,1)]
	y=zz$ycont[c(1:n,1)]
		

	A=sum(x[-(n+1)]*y[-1]-x[-1]*y[-(n+1)])/2
	xb=sum((x[-(n+1)]+x[-1])*(x[-(n+1)]*y[-1]-x[-1]*y[-(n+1)]))/6/A
	yb=sum((y[-(n+1)]+y[-1])*(x[-(n+1)]*y[-1]-x[-1]*y[-(n+1)]))/6/A
	
	c(xb,yb)
			
}


