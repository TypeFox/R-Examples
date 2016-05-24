predictsym.lm <-
function(model,new.sym.data,
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
    for(i in 1:nn) 
      for(j in 1:mm) {
        mins[i,j]<-sym.var(new.sym.data,j)$var.data.vector[i,1]
        maxs[i,j]<-sym.var(new.sym.data,j)$var.data.vector[i,2]
      }
    mins<-as.data.frame(mins)
    colnames(mins)<-new.sym.data$sym.var.names
    maxs<-as.data.frame(maxs)
    colnames(maxs)<-new.sym.data$sym.var.names
    pred.mins<-predict(model,newdata=mins,se.fit=TRUE)
    pred.maxs<-predict(model,newdata=maxs,se.fit=TRUE)    
    Prediction<-data.frame(Minimums=pred.mins$fit, Maximums=pred.maxs$fit)
    return(list(MinPrediction=pred.mins, MaxPredictions=pred.maxs,Fitted=Prediction))
  }
  if(method=='crm') {
    # Center Model
    centers<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(new.sym.data,j)$var.data.vector[i,1]+
                         sym.var(new.sym.data,j)$var.data.vector[i,2])/2    
    centers<-as.data.frame(centers)
    colnames(centers)<-new.sym.data$sym.var.names
    predc<-predict(model$CenterModel,newdata=centers,se.fit=TRUE)
    # Range Model    
    range<-matrix(0,nn,mm)
    for(i in 1:nn) 
      for(j in 1:mm)
        range[i,j]<-(sym.var(new.sym.data,j)$var.data.vector[i,2]-
                       sym.var(new.sym.data,j)$var.data.vector[i,1])/2    
    range<-as.data.frame(range)
    colnames(range)<-new.sym.data$sym.var.names
    predr<-predict(model$RangeModel,newdata=range,se.fit=TRUE)
    res.min<-predc$fit-predr$fit
    res.max<-predc$fit+predr$fit
    Prediction<-data.frame(Minimums=res.min, Maximums=res.max)
    return(list(CenterPrediction=predc,RangePrediction=predr,Fitted=Prediction))    
  }  
}
