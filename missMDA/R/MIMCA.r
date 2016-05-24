MIMCA<-function(X,nboot=100,ncp,coeff.ridge=1,threshold = 1e-06, maxiter = 1000,verbose=FALSE){
  imputeMCA.print<-function (don, ncp, method = c("Regularized", "EM"), row.w = NULL, 
                             coeff.ridge = 1, threshold = 1e-06, seed = NULL, maxiter = 1000,verbose,printm){
    if(verbose){cat(paste(printm,"...",sep=""))}
    res<-imputeMCA(don=don, ncp=ncp , method=method , row.w=row.w, 
                   coeff.ridge=coeff.ridge , threshold = threshold, seed = seed, maxiter = maxiter )
    return(res)
  }
  normtdc<-function(tab.disj,data.na){
    #scale a fuzzy table
    tdc<-tab.disj
    tdc[tdc<0]<-0
    tdc[tdc>1]<-1
    col.suppr<-cumsum(sapply(data.na,function(x){nlevels(x)}))
    tdc<-t(apply(tdc,1,FUN=function(x,col.suppr){
      if(sum(x[1:col.suppr[1]])!=1){x[1:col.suppr[1]]<-x[1:col.suppr[1]]/sum(x[1:col.suppr[1]])}
      for (i in 2:length(col.suppr)){
        x[(col.suppr[i-1]+1):(col.suppr[i])]<-x[(col.suppr[i-1]+1):(col.suppr[i])]/sum(x[(col.suppr[i-1]+1):col.suppr[i]])
      }
      return(x)
    },col.suppr=col.suppr))
    return(tdc)
  }
  draw<-function(tabdisj,Don,Don.na){
    #draw from a scaled fuzzy table
    nbdummy <- rep(1, ncol(Don))
    is.quali <- which(!unlist(lapply(Don, is.numeric)))
    nbdummy[is.quali] <- unlist(lapply(Don[, is.quali, drop = FALSE], 
                                       nlevels))
    vec = c(0, cumsum(nbdummy))
    Donres <- Don
    for (i in is.quali){
      Donres[, i] <- as.factor(levels(Don[, i])[
        apply(tabdisj[,(vec[i] + 1):vec[i + 1]]
              , 1, function(x){
                sample(1:length(x),size=1,prob=x)
              }
        )])
      Donres[, i]<-factor(Donres[, i],levels(Don[,is.quali][,i]))#to keep the order of the levels
    }
    return(don.imp=Donres)
  }
  temp<-if(coeff.ridge==1){"regularized"}else if(coeff.ridge==0){"EM"}else{paste("coeff.ridge=",coeff.ridge)}
  if(verbose){cat("Multiple Imputation using",temp,"MCA using",nboot,"imputed arrays","\n")}
  n<-nrow(X)
  Boot<-matrix(sample(1:n,size=nboot*n,replace=T),n,nboot)
  Weight<-matrix(1/(n*1000),n,nboot,dimnames=list(1:n,paste("nboot=",1:nboot,sep="")))
  Boot.table<-apply(Boot,2,table)
  for(i in 1:nboot){
    Weight[names(Boot.table[[i]]),i]<-Boot.table[[i]]
  }
  Weight<-sweep(Weight,2,STATS=colSums(Weight),FUN="/")
  Weight<-as.data.frame(Weight)
  res.imp<-mapply(Weight,FUN=imputeMCA.print,MoreArgs=list(don=X,ncp=ncp,coeff.ridge= coeff.ridge,method="Regularized",threshold = threshold, maxiter = maxiter,verbose=verbose),printm=as.character(1:nboot),SIMPLIFY=FALSE)
  tdc.imp<-lapply(res.imp,"[[","tab.disj")
  res.comp<-lapply(res.imp,"[[","completeObs")
  tdc.norm<-mapply(FUN=normtdc,tab.disj=tdc.imp,data.na=res.comp,SIMPLIFY=F)
  X.imp<-mapply(FUN=draw,tabdisj=tdc.norm,Don=res.comp,MoreArgs = list(Don.na=X),SIMPLIFY=F)
  res<-list(res.MI=X.imp,
            res.imputeMCA=imputeMCA(X,ncp=ncp,coeff.ridge=coeff.ridge,threshold=threshold,seed=NULL,maxiter=maxiter)$tab.disj,
            call = list(X=X,nboot=nboot,ncp=ncp,coeff.ridge=coeff.ridge,threshold =threshold, seed = NULL, maxiter = maxiter,tab.disj=array(unlist(tdc.imp), dim = c(nrow(tdc.imp[[1]]), ncol(tdc.imp[[1]]), length(tdc.imp))))
  )
  class(res)<-c("MIMCA","list")
  if(verbose){cat("\ndone!\n")}
  return(res)
}