quantile.outliers.trim<-function(x,dfunc=depth.mode,trim=0.25,nb=200,
smo=0.05,ns=0.01,...){
 if (!is.fdata(x)) x=fdata(x)
 dat<-x[["data"]]
 tt<-x[["argvals"]]
 rtt<-x[["rangeval"]]
 n<-nrow(x)
 m<-ncol(x)
 if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
 d=dfunc(x,trim=trim,...)$dep
 rid<-rank(d,ties.method="first")
 num.boot<-floor(trim*n)
 muestra.trim<-x[rid>=num.boot,]
 cuantiles<-numeric(nb)
 vv=var(dat)
 for (i in 1:nb){
        bmuestra<-muestra.trim[sample(1:(n-num.boot),size=n,replace=TRUE),]
       if (smo>0) {bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=n,rep(0,m),vv*smo)}
        d=dfunc(bmuestra,trim=trim,...)$dep
        cuantiles[i]<-quantile(d,probs=ns,type=8)
    }
    return(cuantiles)
}

