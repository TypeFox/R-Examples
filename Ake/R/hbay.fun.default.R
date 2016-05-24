hbay.fun.default <-
function(Vec,x=NULL,...){
y1<-sort(Vec)	
if(is.null(x)){x<-0:max(y1)}	

		vec1=0
		vec2=0
		alp=0.5 
		bet=15  
 		for (i in 1: length(x)){
		    if (x[i]<= y1+1){
			k=seq(0,x[i],by=1)
			vec1[i]=sum ((factorial(y1+1)*(y1^k)*beta(x[i]+alp-k+1,y1+bet-x[i]+1))/(factorial(y1+1-x[i])*factorial(k)*factorial (x[i]-k)*(y1+1)^(y1+1)))
			vec2[i]=sum ((factorial(y1+1)*(y1^k)*beta(x[i]+alp-k,y1+bet-x[i]+1))/(factorial(y1+1-x[i])*factorial(k)*factorial (x[i]-k)*(y1+1)^(y1+1)))
			}
	         else{
                        vec1[i]=0
			vec2[i]=0
			}
                     }
           result=sum(vec1)/sum(vec2)
structure(list(hby=result),class="hbay.fun")
      }
