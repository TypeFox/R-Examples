adaptive.chi.square.main=function(x,B,S,alpha=0.05,prop=0.5,bit=FALSE){
  
  N=length(x)
  M=B
  if (bit==TRUE){
	  say=1
	  i=0
	  biti=seq(0,N,B)
	  M=length(biti)
	  bitis=biti[2:M]
	  bsl=seq(1,N,B)
	  B=matrix(0,(M-1),1)
	  for (i in 1:(M-1)){
  		B[i]=sum(2 ^ (which(as.logical(rev(x[bsl[i]:bitis[i]]))) - 1))
	  }
    x=B
  }    
	
  egitimG=round(N*prop)
  egitimVeri=x[1:egitimG]
  testVeri=x[(egitimG+1):N]

  freq=0
  freq2=0
  
  ss=unique(x)
  if (M>64){
    ss=new('mpfr',unlist(ss))
  }
  for (i in 1:length(ss)){
		freq[i]=sum((ss[i]==egitimVeri)==TRUE)
        
		freq2[i]=sum((ss[i]==testVeri)==TRUE)
   }
   AA=0
   AA2=0
   for (i in 1:S){
		AA[i]=sum(freq==(i-1))
		AA2[i]=sum(freq2==(i-1))	
   }

   statistic=sum((AA2[which(AA2>0)]-AA[which(AA>0)])^2/AA[which(AA>0)])
   pDegeri=1-pchisq(statistic,(S-1))
	
   if (pDegeri<alpha){
		sonucUKK=0
   } else{
		sonucUKK=1
   }
	
   result=list(statistic=statistic,p.value=pDegeri,result.acsq=sonucUKK,name="Achi")
   return(result)
 
}