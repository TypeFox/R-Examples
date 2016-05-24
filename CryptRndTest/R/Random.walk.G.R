Random.walk.G=function(a,B){

  N=length(a)
  b=0
  e=0
  S=0
  b=1-2*a

    e=array(0,floor(N/B)) 
	for (j in 1:floor(N/B)){

		maks=0
		minn=0
		for (i in 1:B){
		  k=1:i
		  S[i]=sum(b[(j-1)*B+k])
   		  if (S[i]>maks){
				maks=S[i]
		  }
		  if (-S[i]>minn){
				minn=-S[i]
		  }
		}
		e[j]=maks+minn
	}
	
	result=list(e=e)
	return(result)
}

