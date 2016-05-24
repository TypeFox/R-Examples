summary.ivmodel <- function(object, ...){
  ivmodel<-object
### print formula
  cat("\nCall:\n", paste(deparse(ivmodel$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("sample size: ", ivmodel$n, "\n", sep="")
  
### first stage regression result
  cat(rep("_", 30), "\n")
  cat("\nFirst Stage Regression Result:\n\n")

  SSM=c(sum(qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj)^2))
  SST=sum(ivmodel$Dadj^2)
  SSE=SST-SSM
  DM=ivmodel$L
  DE=ivmodel$n-ivmodel$p-ivmodel$L
  DT=DM+DE
  
  Fstat=SSM/SSE*DE/DM
  pval=1-pf(Fstat, df1=DM, df2=DE)
  RSquare=SSM/SST
  adjRSquare=1-(1-RSquare)*DT/DE
  RMSE=sqrt(SSE/DE)
  
  cat("F=", Fstat, ", df1=", DM, ", df2=", DE, ", p-value is ", format.pval(pval), "\n", sep="")
  cat("R-squared=", RSquare, ",   Adjusted R-squared=", adjRSquare, "\n", sep="")
  cat("Residual standard error: ", RMSE, " on ", DT, " degrees of freedom\n", sep="")
	  
### Sargan test
  if(ivmodel$L>1){

  cat(rep("_", 30), "\n")
  cat("\nSargan Test Result:\n\n")

  TSLS=sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Yadj))/
       sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj))
  e=ivmodel$Yadj-ivmodel$Dadj*TSLS
  Sargan=sum(qr.fitted(ivmodel$ZadjQR, e)^2)/(sum(e^2)/length(e))
  pval=1-pchisq(Sargan, df=ivmodel$L-1)
  
  cat("Sargan Test Statistics=", Sargan, ", df=", ivmodel$L-1, ", p-value is ", format.pval(pval), "\n", sep="")
  
  }
  
### print TSLS, kClass, LIML, Fuller
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of k-Class Estimators:\n\n")
  printCoefmat(coef(ivmodel), digits = max(3L, getOption("digits") - 3L))
	  
### print AR, ARsens, CLR
  cat(rep("_", 30), "\n")
  cat("\nAlternative tests for the treatment effect under H_0: beta=", ivmodel$beta0, ".\n",sep = "")
  if(!is.null(ivmodel$AR)){
    cat("\nAnderson-Rubin test:\n")
	cat("F=", ivmodel$AR$Fstat, ", df1=", ivmodel$AR$df[1], ", df2=", 
	    ivmodel$AR$df[2], ", p-value=", format.pval(ivmodel$AR$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$AR$ci.info)
  }
  if(!is.null(ivmodel$ARsens)){
    cat("\n\nSensitivity analysis with deltarange [", ivmodel$ARsens$deltarange[1], 
	    ", ", ivmodel$ARsens$deltarange[2], "]:\n")
	cat("non-central F=", ivmodel$ARsens$ncFstat, ", df1=", ivmodel$ARsens$df[1], 
	    ", df2=", ivmodel$ARsens$df[2], ", ncp=", ivmodel$ARsens$ncp, ", p-value=", format.pval(ivmodel$ARsens$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$ARsens$ci.info)
  }
  if(!is.null(ivmodel$CLR)){
    cat("\n\nConditional Likelihood Ratio test:\n")
	cat("Test Stat=", ivmodel$CLR$test.stat, ", p-value=", format.pval(ivmodel$CLR$    p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$CLR$ci.info)
  }
  cat("\n")
}

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