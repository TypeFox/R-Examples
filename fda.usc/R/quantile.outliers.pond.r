quantile.outliers.pond<-function(x,dfunc=depth.mode,nb=200,
smo=0.05,ns=0.01,...){
 if (!is.fdata(x)) x=fdata(x)
 dat<-x[["data"]]
 tt<-x[["argvals"]]
 rtt<-x[["rangeval"]]
 n<-nrow(dat)
 m<-ncol(dat)
 if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
 d=dfunc(x,...)$dep
 cuantiles<-numeric(nb)
 vv=var(dat)
 pr=d/sum(d)
 for (i in 1:nb){
   bmuestra<-x[sample(1:n,size=n,replace=TRUE,prob=pr),]
   if (smo>0) {bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=n,rep(0,m),vv*smo)}
   d=dfunc(bmuestra,...)$dep
   cuantiles[i]<-quantile(d,probs=ns,type=8)
 }
return(cuantiles)
}


