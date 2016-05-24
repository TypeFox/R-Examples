


normalization = function (Xtrain,Xtest=NULL, method = "pqn",ref=NULL) {
    METHODS <- c("pqn", "sum", "median", "sqrt","none")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid normalization method")
    if(method==1) {
   rSXtrain= rowSums(abs(Xtrain), na.rm=T)
   Xtrain = Xtrain / rSXtrain
if (is.null(ref)) ref = apply(Xtrain,2, median, na.rm=T)
 newXtrain =t(t(Xtrain)/ref)
coeXtrain= apply(newXtrain, 1, median, na.rm=T)

if(!is.null(Xtest)){
   rSXtest=rowSums(abs(Xtest), na.rm=T)
Xtest = Xtest / rSXtest
newXtest = t(t(Xtest)/ref)
coeXtest= apply(newXtest,1, median, na.rm=T)
   return(list(newXtrain= Xtrain /coeXtrain, coeXtrain= coeXtrain*rSXtrain,
               newXtest = Xtest / coeXtest, coeXtest= coeXtest*rSXtest)) 

}
return(list(newXtrain= Xtrain /coeXtrain, coeXtrain= coeXtrain*rSXtrain))
}
   
    if(method==2) {
   w1=apply(Xtrain,1,sum)
   data_train=Xtrain/w1
   if(!is.null(Xtest)){
      w2=apply(Xtest,1,function(x) sum(abs(x),na.rm=T))
      data_test=Xtest/w2
      return(list(newXtrain=data_train,coeXtrain=w1,newXtest=data_test,coeXtest=w2))
   }
return(list(newXtrain=data_train,coeXtrain=w1))
}
    
    if(method==3) {
   w1=apply(Xtrain,1,median)
   data_train=Xtrain/w1
   if(!is.null(Xtest)){
      w2=apply(Xtest,1,function(x) median(x,na.rm=T))
      data_test=Xtest/w2
      return(list(newXtrain=data_train,coeXtrain=w1,newXtest=data_test,coeXtest=w2))
   }
return(list(newXtrain=data_train,coeXtrain=w1))
}
 
    if(method==4) {
   w1=apply(Xtrain,1,function(x) sqrt(sum(x^2)))
   data_train=Xtrain/w1
   if(!is.null(Xtest)){
      w2=apply(Xtest,1,function(x) sqrt(sum(x^2,na.rm=T)))
      data_test=Xtest/w2
      return(list(newXtrain=data_train,coeXtrain=w1,newXtest=data_test,coeXtest=w2))
   }
return(list(newXtrain=data_train,coeXtrain=w1))
}
  
    if(method==5) {
   w1=rep(1,nrow(Xtrain))
   data_train=Xtrain
   if(!is.null(Xtest)){
      w2=rep(1,nrow(Xtest))
      data_test=Xtest
      return(list(newXtrain=data_train,coeXtrain=w1,newXtest=data_test,coeXtest=w2))
   }
return(list(newXtrain=data_train,coeXtrain=w1))
}  
}

scaling = function (Xtrain,Xtest=NULL, method = "autoscaling") {
    METHODS <- c("none", "centering", "autoscaling", "rangescaling","paretoscaling")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid scaling method")
    w1=rep(1,nrow(Xtrain))

    if(!is.null(Xtest)){   
    	w2=rep(1,nrow(Xtest))    
        if(method==2) {
            Xtrain <- scale(Xtrain, center=T, scale=F)
            Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=F)
            return(list(newXtrain=Xtrain,newXtest=Xtest))
        }  
        if(method==3) {
        	Xtrain <- scale(Xtrain, center=T, scale=T)
        	Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=attr(Xtrain,"scaled:scale"))	
        }    
        if(method==4) {
        	ran <- apply( Xtrain, 2, function(x) max(x)-min(x) )     
        	Xtrain <- scale(Xtrain, center=T, scale=ran)
        	Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=ran)
        }
        if(method==5) {
        	ssd <- sqrt( apply( Xtrain, 2, sd ) )     
        	Xtrain <- scale(Xtrain, center=T, scale=ssd)
        	Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=ssd)
        }  
        return(list(newXtrain=Xtrain,newXtest=Xtest))  
    }else{
        if(method==2) {Xtrain <- scale(Xtrain, center=T, scale=F) }  
        if(method==3) {Xtrain <- scale(Xtrain, center=T, scale=T)}    
        if(method==4) {
        	ran <- apply( Xtrain, 2, function(x) max(x)-min(x) )     
        	Xtrain <- scale(Xtrain, center=T, scale=ran)
        }
        if(method==5) {
        	ssd <- sqrt( apply( Xtrain, 2, sd ) )     
        	Xtrain <- scale(Xtrain, center=T, scale=ssd)
        }   
        return(list(newXtrain=Xtrain))     	
    }
}







