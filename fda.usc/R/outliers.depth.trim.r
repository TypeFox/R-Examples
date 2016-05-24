
outliers.depth.trim<-function(fdataobj,nb=200,smo=0.05,trim=0.01,quan=0.5,
dfunc=depth.mode,...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")         
 x<-fdataobj[["data"]]
 tt<-fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
 n<-nrow(fdataobj)
 m<-ncol(fdataobj)
# print(dfunc=depth.mode)
 if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
 if (is.null(row.names(fdataobj[["data"]]))) row.names(fdataobj[["data"]])=1:n
    cutoff<-quantile(quantile.outliers.trim(fdataobj,dfunc=dfunc,nb=nb,smo=smo,
    trim=trim,...),probs=quan)
    hay<-1
    outliers<-dep.out<-ite<-c()
    ii<-1
    curvasgood<-fdataobj
    dd<-dfunc(curvasgood,trim=trim,...)     
    modal<-FALSE
    if (!is.null(dd$dist)) {
      modal=TRUE   
      dd<-dfunc(curvasgood,trim=trim,...)
          }
    d<-dd$dep          
    rwn=names(d)=rownames(curvasgood[["data"]])=1:n
    while (hay==1){              
          if (is.null(outliers)){dtotal<-d}
          cutt<-d<cutoff
          fecha<-as.numeric(rownames(curvasgood[["data"]])[cutt])            
          elim<-which(cutt)
          if (length(elim)>0){
             dep.out<-c(dep.out,d[elim])
             curvasgood<-curvasgood[-elim,]
             outliers<-c(outliers,fecha)     
          }    
        if (length(elim)==0 || length(outliers)>n/5){hay<-0}
        else {
            if (modal) {
             mdist<-dd$dist[-elim,-elim]
            class(mdist)<-c("matrix","fdist")        
            dd<-dfunc(curvasgood,trim=trim,metric=mdist,scale=FALSE)
            }
          else dd<-dfunc(curvasgood,trim=trim,...)
          d<-dd$dep 
          }
          ite<-c(ite,rep(ii,length(elim)))
          ii<-ii+1
   }
outliers<-rownames(fdataobj[["data"]])[outliers]    
names(dep.out)<-NULL    
return(list("outliers"=outliers,"dep.out"=dep.out,"iteration"=ite,"quantile"=cutoff,"Dep"=dtotal))
}
