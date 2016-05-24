depth=function(u,x,method="Tukey",approx=FALSE,eps=1e-8,ndir=1000){

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
# Calculates depth of u in sample x.
# u= vector of p real numbers
# x= matrix n by p
# method= which definition of depth to use


  match.arg(method,c("Tukey","Liu","Oja"))
  if(!is.matrix(x)){
    n=length(x)
    if(length(u)>1){stop("Data is univariate, but u is mulativariate.")}
    if(method=="Tukey"){
      return(min(sum(x<=u),sum(x>=u))/n)
    }
    if(method=="Liu"){
      return(sum(outer(u-x,u-x,function(a,b) a*b)<=0)/(n*(n-1)))
    }
    if(method=="Oja"){
      return(.5*(1+mean(abs(u-x)))^(-1))
    }
  }
  p=length(x[1,])
  n=length(x[,1])

  if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.\n")) }
  if(length(u)!=p){ stop("Dimension mismatch between the data and the point u.") }
  if(method!="Tukey"&approx==TRUE){ warning("An approximate algorithm is available only for Tukey's depth. Argument approx=TRUE is ignored.") }
  if(method=="Liu"&p>2){ stop("Liu's depth can be calculated on bivariate datasets only.") }
  if(method=="Tukey"&approx==FALSE&p>3){ warning("Tukey's depth can only be approximated in more than 3 dimensions. Argument approx=F is ignored.") }

# Calculation of Tukey's Median

  if(method=="Tukey"){
  if(p==2){

    ans=.Fortran("fdepth",
      as.numeric(u[1]),
      as.numeric(u[2]),
      as.integer(n),
      as.numeric(x[,1]),
      as.numeric(x[,2]),
      numeric(n),
      integer(n),
      sdep=as.numeric(-1),
      hdep=as.numeric(-1),
	PACKAGE="depth")

      rep=ans$hdep
  }

  if(p==3&&approx==FALSE) {

    zz=.Fortran("stand",
      as.integer(n),
      x=as.numeric(x[,1]),
      y=as.numeric(x[,2]),
      z=as.numeric(x[,3]),
      u=as.numeric(u[1]),
      v=as.numeric(u[2]),
      w=as.numeric(u[3]),
      xn=numeric(n),
      as.numeric(eps),
      err=as.integer(0),
	PACKAGE="depth")

    if(zz$err!=0) {
      if(zz$err>0&&zz$err<10) { stop(message=paste("At least half of the data set share the same value for variable ",zz$err," ."))}
      if(zz$err>10){ stop(message=paste("Variable ",zz$err-10," has null covariance."))}
    }

    ans=.Fortran("depth3",
      as.integer(n),
      as.numeric(zz$u),
      as.numeric(zz$v),
      as.numeric(zz$w),
      as.numeric(zz$x),
      as.numeric(zz$y),
      as.numeric(zz$z),
      numeric(n),
      integer(n),
      numeric(n),
      numeric(n),
      as.numeric(eps),
      ndim=as.integer(0),
      ndep=as.integer(0),
	PACKAGE="depth")

    rep=ans$ndep/n
  }

  if(approx==TRUE||p>3){

    zz=.Fortran("standpd",
      as.integer(n),
      as.integer(p),
      n=as.integer(n),
      p=as.integer(p),
      x=as.numeric(x),
      u=as.numeric(u),
      as.numeric(eps),
      err=integer(p),
      ndep=as.numeric(n),
	PACKAGE="depth")

    if(zz$ndep==0||zz$p==0) {
      cat("\nSample is singular.\n")
      stop()
    }


    z=.Fortran("hdepth",
      as.integer(n),
      as.integer(p),
      nnp=as.integer(p),
      as.integer(ndir),
      as.integer(n),
      as.integer(p),
      as.numeric(zz$x),
      integer(p),
      as.numeric(zz$u),
      numeric(p),
      numeric(p*p),
      numeric(p),
      numeric(p*p),
      numeric(p),
      as.numeric(eps),
      err=as.numeric(0),
      err2=numeric(ndir),
      ndep=as.integer(0),
      as.integer(0),
	PACKAGE="depth")

    if (z$err==-1) { warning(message=paste("\nThe dataset is singular.\nIts dimension was reduced to ",z$nnp, " dimensions.")) }
    if (z$err==-2) { stop("Error in dimension reduction: eigenvectors are not independent.") }
    for (i in 1:ndir){
      if (z$err2[i]==-1) { stop(message=paste("\nNo eigenvalue for datum ",z$err2[i])) }
      if (z$err2[i]==-2) { stop(message=paste("\nNull eigenvector for datum ",z$err2[i])) }
      if (z$err2[i]>0) {   stop(message=paste("\nError ",z$err2[i], " in the calculation of eigenvalues.")) }
    }
    
    rep=z$ndep/n
    }
  }

# Calculation of Liu's median
  
  if(method=="Liu"){
    ans=.Fortran("fdepth",
      as.numeric(u[1]),
      as.numeric(u[2]),
      as.integer(n),
      as.numeric(x[,1]),
      as.numeric(x[,2]),
      numeric(n),
      integer(n),
      sdep=as.numeric(-1),
      hdep=as.numeric(-1),
	PACKAGE="depth")

      rep=ans$sdep
  
  }
  
  if(method=="Oja"){
    ans=.Fortran("ojadepth",
      as.numeric(x),
      as.numeric(u),
      as.integer(p),
      as.integer(n),
      depth=as.numeric(-1),
	PACKAGE="depth")
      
    rep=ans$depth
  
  }
  
  rep
}
