sym.interval.pca <-
function(sym.data,method=c('classic','tops','centers')) {
  idn<-all(sym.data$sym.var.types==sym.data$sym.var.types[1])
  if(idn==FALSE) 
    stop("All variables have to be of the same type")     
  method<-match.arg(method)  
  if(method=='classic') {
    if((sym.data$sym.var.types[1]!='$C')&&(sym.data$sym.var.types[1]!='$I')) 
      stop("Variables have to be continuos or Interval")    
    if(sym.data$sym.var.types[1]=='$C')
       res<-PCA(sym.data$data, scale.unit=TRUE, ncp=sym.data$M, graph = FALSE)
    else
      if(sym.data$sym.var.types[1]=='$I') {
        nn<-sym.data$N
        mm<-sym.data$M
        centers<-matrix(0,nn,mm)
        centers<-as.data.frame(centers)
        rownames(centers)<-sym.data$sym.obj.names
        colnames(centers)<-sym.data$sym.var.names
        for(i in 1:nn) 
          for(j in 1:mm)
            centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                             sym.var(sym.data,j)$var.data.vector[i,2])/2
        res<-PCA(centers, scale.unit=TRUE, ncp=sym.data$M, graph = FALSE)
      }
    return(res)
  } 
  if(method=='centers') {
    nn<-sym.data$N
    mm<-sym.data$M
    centers<-matrix(0,nn,mm)
    centers.stan<-matrix(0,nn,mm)
    min.stan<-matrix(0,nn,mm)
    max.stan<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                         sym.var(sym.data,j)$var.data.vector[i,2])/2
    # Standarized
    for(i in 1:nn) 
      for(j in 1:mm)
        centers.stan[i,j]<-(centers[i,j]-mean(centers[,j]))/
                           (sd(centers[,j])*sqrt((nn-1)/nn))    
    # Min-Max
    for(i in 1:nn){
      for(j in 1:mm) {
        min.stan[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]-
                          mean(centers[,j]))/(sd(centers[,j])*sqrt((nn-1)/nn))    
        max.stan[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,2]-
                      mean(centers[,j]))/(sd(centers[,j])*sqrt((nn-1)/nn))
      }
    }    
    # Correlation Centers Matrix
    R<-t(centers.stan)%*%centers.stan    
    svd<-eigen(R)
    sym.comp<-sym.data
    # Interval Principal Components
    for(i in 1:nn){
      posd<-1
      for(j in 1:mm) {
        smin<-0
        smax<-0
        for(k in 1:mm) {
          if(svd$vectors[k,j]<0)
            smin<-smin+max.stan[i,k]*svd$vectors[k,j]
          else
            smin<-smin+min.stan[i,k]*svd$vectors[k,j]
          if(svd$vectors[k,j]<0)
            smax<-smax+min.stan[i,k]*svd$vectors[k,j]
          else
            smax<-smax+max.stan[i,k]*svd$vectors[k,j]
        }
        sym.comp$meta[i,sym.comp$sym.var.starts[j]]<-smin
        sym.comp$meta[i,sym.comp$sym.var.starts[j]+1]<-smax
        sym.comp$data[i,posd]<-smin
        sym.comp$data[i,posd+1]<-smax        
        posd<-posd+2
      }
    }
    pos<-1
    for(j in 1:mm) {
      comp.name<-paste("C", j, sep = "")
      sym.comp$sym.var.names[j]<-comp.name
      comp.name<-paste("Min.C", j, sep = "")
      colnames(sym.comp$data)[pos]<-comp.name
      comp.name<-paste("Max.C", j, sep = "")
      colnames(sym.comp$data)[pos+1]<-comp.name
      pos<-pos+2
      comp.name<-paste("Min.C", j, sep = "")
      colnames(sym.comp$meta)[sym.comp$sym.var.starts[j]]<-comp.name
      comp.name<-paste("Max.C", j, sep = "")
      colnames(sym.comp$meta)[sym.comp$sym.var.starts[j]+1]<-comp.name
    }
    # Interval Principal Correlations
    svdV<-matrix(0,nn,nn)
    for(i in 1:nn) {
      for(j in 1:mm) {
        ss<-0
        for(k in 1:mm) {
          ss<-ss+centers.stan[i,k]*svd$vectors[k,j]
        }        
        svdV[i,j]<-(1/sqrt(svd$values[j]))*ss
      }
    }    
    IPrinCorre<-matrix(0,mm,2*mm)
    for(i in 1:mm){
      pcol<-1
      for(j in 1:mm) {
        smin<-0
        smax<-0
        for(k in 1:nn) {
          if(svdV[k,j]<0)
            smin<-smin+(1/sqrt(nn))*max.stan[k,i]*svdV[k,j]
          else
            smin<-smin+(1/sqrt(nn))*min.stan[k,i]*svdV[k,j]
          if(svdV[k,j]<0)
            smax<-smax+(1/sqrt(nn))*min.stan[k,i]*svdV[k,j]
          else
            smax<-smax+(1/sqrt(nn))*max.stan[k,i]*svdV[k,j]
        }
        IPrinCorre[i,pcol]<-smin
        IPrinCorre[i,pcol+1]<-smax
        pcol<-pcol+2
      }
    }
    IPrinCorre<-as.data.frame(IPrinCorre)
    rownames(IPrinCorre)<-sym.data$sym.var.names
    return(list(Sym.Components=sym.comp,Sym.Prin.Correlations=IPrinCorre))
  }  
  return(TRUE)
}
