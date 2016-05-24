TBT.criticalValue=function(m, k, alpha=0.01, cdf=FALSE, exact=TRUE){
# If exact is TRUE it will use gmp package to calculate Stirling number of the second kind. 
# However, this takes too long time!
  j=1
  deneme=Inf
  while ((deneme==Inf)|(deneme==0)){
    j=j+1
    deneme=exp(sum(log(((2^m)-j+1):(2^m)))+Strlng2(k,j)-k*log(2^m))
  }
  bsl=(j-1)
  topl=deneme
  prob=0
  value=0
  critical.value=0
  if (cdf==TRUE){
    say=0
    for (j in bsl:(k-1)){
      hesap=exp(sum(log(((2^m)-j+1):(2^m)))+Strlng2(k,j,log=TRUE)-k*log(2^m))*(2*m+log(k)) # here a correction factor is applied
      if (exact==TRUE){
        while (is.nan(hesap)==TRUE){ #if strlng2  cannot be calculated do the following
          hesap=exp(sum(log(((2^m)-j+1):(2^m)))+gmp::Stirling2(k,j)-k*log(2^m))
        }
      }
      topl=topl+hesap  
      if ((topl>1)|(is.nan(topl)==TRUE)){
        topl=1
      }
      say=say+1
      prob[say]=topl
      value[say]=j
      if (topl>(1-alpha)){
        critical.value=j  
      }
    }
    prob[say]=1
    value[say]=k
    result=list(prob=prob,value=value,critical.value=critical.value)
  }else{ #only finds the critical value
    devam2=1
    j=bsl
    while ((devam2==1)&(j<=k)){    
      hesap=exp(sum(log(((2^m)-j+1):(2^m)))+Strlng2(k,j,log=TRUE)-k*log(2^m))*(2*m+log(k)) # here a correction factor is applied
      if (exact==TRUE){
        while (is.nan(hesap)==TRUE){ #if strlng2  cannot be calculated do the following
          hesap=exp(sum(log(((2^m)-j+1):(2^m)))+gmp::Stirling2(k,j)-k*log(2^m))
        }
      }
      topl=topl+hesap   
      if ((topl>1)|(is.nan(topl)==TRUE)){
        topl=1
      }
      if (topl>(1-alpha)){
        critical.value=j  
        devam2=0
      }
      j=j+1
    }
    result=list(critical.value=critical.value)
  }  
  
  return(result)
}