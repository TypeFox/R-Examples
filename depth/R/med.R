med=function(x,method="Tukey",approx=FALSE,eps=1e-8,maxit=200,mustdith=FALSE,maxdith=50,dithfactor=10,factor=0.8,nstp=NULL,ntry=NULL,nalt=NULL,ndir=1000){

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

  match.arg(method,c("Tukey","Liu","Oja","Spatial","CWmed"))
  if(is.null(dim(x))){ m=median(x)
    return(list(median=m,depth=depth(m,x,method=method))) 
  }
  p=length(x[1,])
  n=length(x[,1])
  if(p==1){ return(median(x)) }
  if(is.null(nstp)) nstp=floor(5*n^0.3*p)
  if(is.null(ntry)) ntry=10*(p+1)
  if(is.null(nalt)) nalt=4*(p+1)
  
  if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.\n")) }
  if(p!=2 & method=="Liu"){ stop("Liu's median can be calculated on bivaraite data sets only") }
  if(p!=2 & method=="Oja"){ stop("Oja's median can be calculated on bivaraite data sets only") } 
  if(p>2 & method=="Tukey" & approx==FALSE){ 
    warning("Tukey's median can be calculated exactly on bivariate samples only.") 
    approx=TRUE
  }
  
  if(method=="Tukey"){
	if (p==2&&approx==FALSE) {

		maxnum=floor(4*n*sqrt(n)+1)

		zz=.Fortran("halfmed",
			as.numeric(x[,1]),
			as.numeric(x[,2]),
			as.integer(n),
			dpth=integer(1),
			med=numeric(2),
			numeric(n*(n-1)/2),
			numeric(n*(n-1)/2),
			integer(1),
			numeric(1),
			integer(1),
			as.integer(0),
			as.integer(maxnum),
			err=integer(1),
			as.numeric(eps),
			as.numeric(dithfactor),
			as.integer(maxdith),
			as.integer(mustdith),
			missing=integer(1),
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
	if (zz$err==-4) { warning("A non critical numerical error occurred during the calculation of a depth contour.") }
	if (zz$err>0) { warning(message=paste("CRITICAL numerical error in calculating contour",zz$err,". The median is probably wrong.")) }

	rep=list(median=zz$med,depth=zz$dpth/n)
	
	}

	if(p>2||approx==TRUE) {
	
		zz=.Fortran("deeplocstand",
			as.integer(n),
			as.integer(p),
			n=as.integer(n),
			p=as.integer(p),
			x=as.numeric(as.matrix(x)),
			xn=numeric(n),
			as.numeric(eps),
			locsca=numeric(p*2),
			err=integer(p),
			PACKAGE="depth")

			for(i in 1:p) {
				if(zz$err[i]==-1){ warning("Variable ",i," have null variance. The number of dimensions was reduced.") }
				if(zz$err[i]==-2){ stop(message=paste("At least half of the data set share the same value for variable ",zz$err," ."))}
			}

		z=.Fortran("deepest",
			as.integer(n),
			as.integer(p),
			as.integer(zz$n),
			as.integer(zz$p),
			as.integer(ndir),
			as.numeric(zz$x),
			as.numeric(eps),
			nddpst=integer(1),
			dpstM=numeric(p),
			as.integer(nstp),
			as.integer(ntry),
			as.integer(nalt),
			err=integer(floor(ndir/4)),
			errc=integer(1),
			as.integer(floor(ndir/4)),
			PACKAGE="depth")
		if(z$errc>0){  warning(message=paste("No more improvement after nstep iterations; ntry = ",z$errc)) }
		if(z$errc<0){  warning(message=paste("Reach maximum number of iterations: nstep = ",-z$errc)) }
			for(i in 1:floor(ndir/4)) {
				if(z$err[i]>0) { warning(message=paste("Error ",z$err[i]," in the calculation of eigenvectors."))}
				if(z$err[i]==-1){ warning(message=paste("No null eigenvalue for sample ",z$err[i]) )}
				if(z$err[i]==-2){ warning(message=paste("Null eigenvector for sample ",z$err[i]))}
			}
		

		locsca=matrix(zz$locsca,p,2)

		for (i in 1:p) {
			z$dpstM[i]=z$dpstM[i]*locsca[i,2]+locsca[i,1]
		}
		rep=list(median=z$dpstM,depth=z$nddpst/n)
	}
  
  }
  
  if(method=="Liu"){
	ans=.Fortran("liumed",
			as.numeric(x[,1]),
			as.numeric(x[,2]),
			as.integer(n),
			dpth=as.integer(0),
			med=numeric(2),
			PACKAGE="depth")
	rep=list(median=ans$med,depth=ans$dpth/(n*(n-1)*(n-2)/6))  
  }
  
  if(method=="Oja"){

	ans=.Fortran("ojamed",
			as.numeric(x[,1]),
			as.numeric(x[,2]),
			as.integer(n),
			integer(n),
			xmed=numeric(1),
			ymed=numeric(1),
			err=integer(1),
			as.numeric(eps),
			PACKAGE="depth")
		if(ans$err==-1) { stop("Colinearity detected in the data.") }
		if(ans$err>0) {	stop(paste("Fatal error : ifault = ", ans$err,". (See AS277).")) }
	omed=c(ans$xmed,ans$ymed)
	rep=list(median=omed,depth=depth(omed,x,method="Oja"))
  
  }
  
  if(method=="Spatial"){

  # Algorithm AS143 is faster; if it fails, algorithm AS78 is used instead.


	ans=.Fortran("medctr",
			as.integer(p),
			as.integer(n),
			as.numeric(x),
			numeric(1),
			numeric(1),
			med=numeric(p),
			err=integer(1),
			as.numeric(eps),
			it=as.integer(maxit),
			PACKAGE="depth")

			if(ans$err==-2||ans$err==1) { 

			ans=.Fortran("medctr78",
				as.numeric(x),
				med=numeric(p),
				as.integer(n),
				as.integer(p),
				integer(1),
				err=integer(1),
				PACKAGE="depth")
			}
	rep=list(median=ans$med)  
  }
  
  if(method=="CWmed"){
    rep=list(median=apply(x,2,median))
  }
  
  rep
}
