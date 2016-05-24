confint.ivmodel <- function(object, parm,level=NULL, ...){
  ivmodel = object
  alpha=ivmodel$alpha
  if(!is.null(level)){
    if(!is.numeric(level) | level>1 | level<0){
	  print("Wrong input of confidence level!")
	  return()
	}else{
	  alpha=1-level
	}
  }
  temp<-coef(ivmodel)
  result<-matrix(NA, ncol=2, nrow=nrow(temp))
  colnames(result)<-c(paste(round(alpha/2*100, digits=2), "%", sep=""),
                      paste(round((1-alpha/2)*100, digits=2), "%", sep=""))
  rownames(result)<-rownames(temp)
  index<-rownames(temp)=="k-class"
  rownames(result)[index]<-paste("k-class", round(temp[index, 1], digits=2))
  result[,1]<-temp[,2]+temp[,3]*qt(alpha/2, df=ivmodel$n-ivmodel$p-1)
  result[,2]<-temp[,2]+temp[,3]*qt(1-alpha/2, df=ivmodel$n-ivmodel$p-1)

  if(!is.null(ivmodel$AR)){
    if(alpha!=ivmodel$alpha){
	  temp=AR.test(ivmodel, beta0=ivmodel$beta0, alpha=alpha)
	}else{
	  temp=ivmodel$AR
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"AR"
	}else{
	  rownames(temp$ci)<-paste("AR(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  if(!is.null(ivmodel$ARsens)){
    if(alpha!=ivmodel$alpha){
	  temp=ARsens.test(ivmodel, beta0=ivmodel$beta0, alpha=alpha, 
	                   deltarange=ivmodel$ARsens$deltarange)
	}else{
	  temp=ivmodel$ARsens
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"Sensitivity Analysis"
	}else{
	  rownames(temp$ci)<-paste("Sensitivity Analysis(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  if(!is.null(ivmodel$CLR)){
    if(alpha!=ivmodel$alpha){
	  temp=CLR(ivmodel, beta0=ivmodel$beta0, alpha=alpha)
	}else{
	  temp=ivmodel$CLR
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"CLR"
	}else{
	  rownames(temp$ci)<-paste("CLR(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  return(result) 
}