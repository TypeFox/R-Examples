Estimator=function(output,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){

  if(!is.list(output)){stop("output must be a list.")}
  if(any(is.na(output))){stop("There are missing values in output.")}
  
  if(!is.vector(pi)){stop("pi must be a vector.")}
  if(any(is.na(pi))){stop("There are missing values in pi.")}
  if(any((pi<=0)|(pi>1))){stop("There are invalid values in pi.")}         
  n=length(pi)
  if(n!=length(output$TransformedVariable)){stop("The lengths of pi and transformed variable are different.")}
  if(n!=length(output$TransformedVariance)){stop("The lengths of pi and transformed variance are different.")}
  
  if(!is.character(type)){stop("type must be a character")}
  if((type!="total")&(type!="mean")){stop("The value of type must be total or mean.")}  
  
  if((cl<=0)|(cl>=1)){stop("The value of cl must be in the range (0,1).")}
   
  zalpha=qnorm(1-(1-cl)/2)
  
  if(type=="total"){
    e=sum(output$TransformedVariable/pi) 
    ve=sum(output$TransformedVariance/pi)+varest(output$TransformedVariable,pik=pi)
    ci=c((e-zalpha*sqrt(ve)),(e+zalpha*sqrt(ve)))
  }
  if(type=="mean"){
    if((is.null(N))&(is.null(pij))){
      N=round(sum(1/pi))
      e=(1/N)*sum(output$TransformedVariable/pi)
      warning("To calculate the estimated variance is needed or the size of the population or the second-order inclusion probabilities matrix.")
      out=list(Estimation=e)
      return(out)
    }
    if(!is.null(N)){
      if(length(N)!=1){stop("N must be a scalar.")}   
      if(N<0){stop("N must be a positive number.")}
    
      ve=(1/N^2)*(sum(output$TransformedVariance/pi)+varest(output$TransformedVariable,pik=pi))
    }else{
      if(!is.matrix(pij)){stop("pij must be a matrix.")}     
      if(nrow(pij)!=ncol(pij)){stop("pij is not a square matrix.")} 
      if(any(is.na(pij))){stop("There are missing values in pij.")}
      if(any((pij<=0)|(pij>1))){stop("There are invalid values in pij")} 
      if(n!=nrow(pij)){stop("The lengths of pi and pij are different.")} 
      
      N=round(sum(1/pi))
      ve=(1/N^2)*(sum(output$TransformedVariance/pi)+vartaylor_ratio(output$TransformedVariable,rep(1,n),pij)$estvar)
    }
    e=(1/N)*sum(output$TransformedVariable/pi)
    ci=c((e-zalpha*sqrt(ve)),(e+zalpha*sqrt(ve)))
  }
  
  if(ve<0){warning("The variance estimation can not be negative.")}
  
  out=list(Estimation=e,Variance=ve,ConfidenceInterval=ci)
  return(out)
}