sym.glm <-
function(sym.data,response=1,method=c('cm','crm'),alpha=1,
                  nfolds=10,grouped=TRUE) {
  idn<-all(sym.data$sym.var.types==sym.data$sym.var.types[1])
  if(idn==FALSE) 
    stop("All variables have to be of the same type")    
  method<-match.arg(method)
  nn<-sym.data$N
  mm<-sym.data$M  
  if(method=='cm') {
    centers<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                         sym.var(sym.data,j)$var.data.vector[i,2])/2    
    model<-cv.glmnet(centers[,-response],centers[,response],nfolds=nfolds,grouped=grouped,
                     alpha=alpha) 
    return(model)
  }  
  if(method=='crm') {
    ## Center Model
    centers<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                         sym.var(sym.data,j)$var.data.vector[i,2])/2    
    modelc<-cv.glmnet(centers[,-response],centers[,response],nfolds=nfolds,grouped=grouped,
                     alpha=alpha) 
    # Range Model    
    range<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        range[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,2]-
                       sym.var(sym.data,j)$var.data.vector[i,1])/2 
    modelr<-cv.glmnet(range[,-response],range[,response],nfolds=nfolds,grouped=grouped,
                      alpha=alpha) 
    return(list(CenterModel=modelc,RangeModel=modelr))
  }  
}
