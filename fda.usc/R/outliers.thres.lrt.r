outliers.thres.lrt<-function(fdataobj,nb=200,smo=0.05,trim=0.10,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
names<-fdataobj[["names"]]
rtt<-fdataobj[["rangeval"]]
    n<-nrow(fdataobj)
    m<-ncol(fdataobj)
    if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
    maximos<-c()
    aux<-c()
    remuestras.boot<-matrix(NA,nrow=n,ncol=m)
    pb=txtProgressBar(min=0,max=nb,style=3)
    for (i in 1:nb){
        setTxtProgressBar(pb,i-0.5)
        bmuestra<-fdataobj[sample(1:n,size=n,replace=TRUE),]
        if (smo>0) {bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=n,rep(0,m),var(fdataobj[["data"]])*smo)}
        auxmean<-func.trim.mode(bmuestra,trim=trim,...)
        auxdt<-sqrt(as.vector(func.trimvar.mode(bmuestra,trim=trim,...)[["data"]]))
        d<-matrix(NA,nrow=n,ncol=m)
        for (j in 1:m){d[,j]<-1-abs(.5-rank(bmuestra[,j][["data"]],ties.method="average")/n)}
#        ans<-apply(d,1,sum)
        ans<-rowSums(d)
        rid<-rank(ans,ties.method="first")
        bmuestra.trim<-bmuestra[rid>=floor(trim*n),]
        for (j in 1:(n-floor(trim*n)))
        {aux[j]<-metric.lp(bmuestra.trim[j,][["data"]]/auxdt,auxmean[["data"]]/auxdt,...)}
        maximos[i]<-as.numeric(max(aux))
        setTxtProgressBar(pb,i)
    }
    close(pb)
    ans<-as.numeric(quantile(maximos,.99))
    ans
}



