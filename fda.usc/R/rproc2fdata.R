rproc2fdata=function(n,t=NULL,mu=rep(0,length(t)),sigma=1,
                     par.list=list("scale"=1,"theta"=.2*diff(rtt),"H"=.5),
                     norm=FALSE,verbose=FALSE,...) {
sigma2<-sigma
if (is.fdata(mu)){   
  if (!is.null(t) & verbose) warnings("The argvals of argument t are ignored, the function uses the argvals(mu)")
  t<-mu$argvals
  p=length(t)   
  rtt<-mu$rangeval
  mu<-drop(mu$data[1,])
 } 
else {
  if (is.null(t)) stop("At least t or fdata class mu must be provided")
  p=length(t)
  if (p!= length(mu)) stop("t and mu must have the same length")
  rtt=range(t)
  }
if (is.null(par.list$mu0)) par.list$mu0<-0
if (is.null(par.list$theta)) par.list$theta<-1/(3*(diff(rtt)))
#if (is.null(par.list$theta)) par.list$theta<-1/3
if (is.null(par.list$scale)) par.list$scale<-1
if (is.character(sigma)) {
 type.proc<-c("brownian","wiener","OU","OrnsteinUhlenbeck","vexponential","fbrownian")
 if (!is.element(sigma,type.proc)) stop("Error in sigma argument label")
 ss<-which(sigma==type.proc)
   sigma2<-sigma
  sigma=switch(sigma,"brownian"=par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "wiener"= par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "OU"= {#m<-100;t2<-t+m
         par.list$scale/(2*par.list$theta)*outer(t,t,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
"OrnsteinUhlenbeck"={#m<-100;#t2<-t+m
            par.list$scale/(2*par.list$theta)*outer(t,t,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
             "vexponential"=par.list$scale*outer(t,t,function(u,v){
             exp(-abs(u-v)/par.list$theta)}),
			"fbrownian"=0.5*abs(par.list$scale)^par.list$H*outer(t,t,function(u,v){abs(u)^(2*par.list$H)+abs(v)^(2*par.list$H)-abs(u-v)^(2*par.list$H)})
             )
  sigma<-t(sigma)

}
else {
 sigma2<-"Gaussian"
 if   (is.matrix(sigma)) {
  if (dim(sigma)[2]!=p) stop("Error in sigma argument")
 }
 else if (length(sigma)==1 | length(sigma)==p) sigma<-diag(p)*sigma
 else stop("Error in sigma argument")
}
C=svd(t(sigma))
L=C$u%*%diag(sqrt(C$d))
X=matrix(rnorm(n*p),ncol=p)
X=t(L%*%t(X))
X=sweep(X,2,mu,"+")
X=fdata(X,t,rtt,names=list("main"=paste(sigma2," process",sep="")))
if (norm) {
if (sigma2[1]=="brownian") print("The normalization is not done")
else{
        no <- norm.fdata(X,...)
        X$data <- sweep(X$data, 1, drop(no), "/")
        }
}
return(X)
}

#{#rgenfdata
rcombfdata=function(n = 10, fdataobj, mu,
                    sdarg = rep(1,nrow(fdataobj)),
                    norm = 1){
  if (class(fdataobj)!="fdata") 
    stop("Argument fdataobj must be of class fdata")
  tt <- argvals(fdataobj)
  if (missing(mu)) 
    mu=fdata(rep(0,ncol(fdataobj)),argvals=tt)
  nr=nrow(fdataobj)
  xx=matrix(rnorm(n*nr),ncol=nr)
  xx=sweep(xx,1,apply(xx,1,function(v){sqrt(sum(v^2))})/norm,"/")
  xx=sweep(xx,2,sdarg,"*")
  res=xx%*%fdataobj[["data"]]
  res=fdata(sweep(res,2,mu[["data"]],"+"),argvals=tt)
  return(res)
}
#}

#{#gridfdata
gridfdata=function(coef,fdataobj,mu){
  if (class(fdataobj)!="fdata") stop("Argument fdataobj must be of class fdata")
  nr=nrow(fdataobj)
  tt <- argvals(fdataobj)
  if (missing(mu)) 
    mu=fdata(rep(0,ncol(fdataobj)),argvals=tt)  
  coef=as.matrix(coef)
  if (ncol(coef)!=nr) stop("Argument coef must be a matrix with ncol(coef)==nrow(fdataobj)")
  res=coef%*%fdataobj[["data"]]
  res=fdata(sweep(res,2,mu[["data"]],"+"),argvals=tt)
  return(res)
}
#}

