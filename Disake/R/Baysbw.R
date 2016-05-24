Baysbw <-
function (Vec){
###########################################################################################################
#    The bayesian approach is used only with the binomial kernel.
#==========================================================================================================                                     
# INPUT:
#   "Vec"	: sample of data   
# OUTPUT:
#    Returns the bandwidth computed using the local Bayesian approach.
########################################################################################################### 
y1<-sort(Vec)		
x<-0:max(y1)
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
           return(sum(vec1)/sum(vec2))
      }
