predictsym.glm <-
function(model,new.sym.data,response=1,
              method=c('cm','crm')) {
  idn<-all(new.sym.data$sym.var.types==new.sym.data$sym.var.types[1])
  if(idn==FALSE) 
    stop("All variables have to be of the same type")    
  method<-match.arg(method)
  nn<-new.sym.data$N
  mm<-new.sym.data$M  
  if(method=='cm') {
    mins<-matrix(0,nn,mm)
    maxs<-matrix(0,nn,mm)
    for(i in 1:nn) {
      for(j in 1:mm) {
        mins[i,j]<-sym.var(new.sym.data,j)$var.data.vector[i,1]
        maxs[i,j]<-sym.var(new.sym.data,j)$var.data.vector[i,2]
      }
    }  
    pred.mins<-predict(model,newx=mins[,-response], s="lambda.min")   
    pred.maxs<-predict(model,newx=maxs[,-response], s="lambda.min") 
    Prediction<-data.frame(Minimums=pred.mins, Maximums=pred.maxs)
    return(Prediction)
  }  
  if(method=='crm') {
    # Center Model
    centers<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(new.sym.data,j)$var.data.vector[i,1]+
                         sym.var(new.sym.data,j)$var.data.vector[i,2])/2    
    predc<-predict(model$CenterModel,newx=centers[,-response], s="lambda.min")   
    # Range Model    
    range<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        range[i,j]<-(sym.var(new.sym.data,j)$var.data.vector[i,2]-
                       sym.var(new.sym.data,j)$var.data.vector[i,1])/2    
    predr<-predict(model$RangeModel,newx=range[,-response], s="lambda.min") 
    res.min<-predc-predr
    res.max<-predc+predr
    Prediction<-data.frame(Minimums=res.min, Maximums=res.max)
    return(Prediction)    
  }  
}
