book.stack.main=function(x,B,k=2,alpha=0.05,bit=FALSE){
  
  if (bit==TRUE){
	  N=length(x)	
	  say=1
	  i=0
  	biti=seq(0,N,B)
	  M=length(biti)
  	bitis=biti[2:M]
	  bsl=seq(1,N,B)
	  A=matrix(0,(M-1),1)
	  for (i in 1:(M-1)){
		  A[i]=sum(2 ^ (which(as.logical(rev(x[bsl[i]:bitis[i]]))) - 1))
    }
    x=A
  }
  
  alfabe=unique(x)

  M=length(alfabe)
  N=length(x)
  e=array(0,k)
  nu=matrix(1:(M+1),nrow=(N+1), ncol=(M+1))  
  
  #nu[1,]=1:(M+1)  
  m=ceiling(max(x)/k)
  for (i in 2:(N+1)){
      if (i%%1000==0){
        print("burada")
      }
      nu[i,1]=x[i-1]
      a=ceiling(which(nu[i-1,]==x[i-1])/m)
      if (length(a)!=0){
          e[a]=e[a]+1
          nu[i,2:(M+1)]=nu[i-1,which(nu[i-1,]!=x[i-1])]          
      } else {
          a.c=ceiling(x[i-1]/m)
          e[a.c]=e[a.c]+1
          gecis=nu[i-1,which(nu[i-1,]!=x[i-1])]
          nu[i,2:(M+1)]=gecis[1:(length(gecis)-1)]          
      }
  }
	
  kiKare=sum((e-m)^2/m)
  pDegeri=1-pchisq(kiKare,(k-1))
 
  if (pDegeri<alpha){
	sonucKR=0
  } else{
	sonucKR=1
  }
  
  result=list(statistic=kiKare,p.value=pDegeri,BS.result=sonucKR,name="BS")
  return(result)
}