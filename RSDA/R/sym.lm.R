sym.lm <-
sym.lm<-function(formula,sym.data,method=c('cm','crm')) {
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
      centers<-as.data.frame(centers)
      colnames(centers)<-sym.data$sym.var.names
      model<-lm(formula,data=centers)
      return(model)
    }
    if(method=='crm') {
      # Center Model
      centers<-matrix(0,nn,mm)
      for(i in 1:nn) 
        for(j in 1:mm)
          centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                           sym.var(sym.data,j)$var.data.vector[i,2])/2    
      centers<-as.data.frame(centers)
      colnames(centers)<-sym.data$sym.var.names
      modelc<-lm(formula,data=centers)
      # Range Model    
      range<-matrix(0,nn,mm)
      for(i in 1:nn) 
        for(j in 1:mm)
          range[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,2]-
                         sym.var(sym.data,j)$var.data.vector[i,1])/2    
      range<-as.data.frame(range)
      colnames(range)<-sym.data$sym.var.names
      modelr<-lm(formula,data=range)
      return(list(CenterModel=modelc,RangeModel=modelr))
    }  
  }