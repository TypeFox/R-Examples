Random.walk.Y=function(a,B){
  N=length(a)
  b=1-2*a

  e=array(0,floor(N/B)) 
	for (j in 1:floor(N/B)){
	  S=0
   	  for (i in 1:B){
		  k=1:i
		  S[i]=sum(b[(j-1)*B+k])
   		  if (abs(S[i])>e[j]){
				e[j]=abs(S[i])
		  }
 		}
	}
	result=list(e=e)
	return(result)
}

