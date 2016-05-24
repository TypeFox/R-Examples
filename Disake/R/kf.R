kf <-
function(x,t,h,ker,a=1,c=2){
###########################################################################################################
# INPUTS:
#   "x" 	: The target.
#   "t" 	: The single value where the function is computed.
#   "h" 	: The bandwidth. It is in (0,1]for the binomial kernel.
#   "ker" 	: The kernel:  "dirDU" DiracDU,"bino" Binomial,"triang" Discrete Triangular.	 
#   "a" 	: The arm is used only for the discrete triangular distribution.
#   "c" 	: The number of categories in the Aitchison and Aitken kernel and is used only for DiracDU
# OUTPUT:
# Returns the discrete associated kernel estimation at t.
###########################################################################################################
 dtrg=function(x,t,h,a){

if (a==0)
	{
			 result <- t
       			 Logic1 <- (t==x) 
			 Logic0 <- (t!=x) 
       			 result[Logic1]=1
        		 result[Logic0]=0					
        		return(result)
			}
  
else

      {	
		
			 u=0:a;
 			 u=sum(u^h)			 
			 D=(2*a+1)*(a+1)^h -2*u                 
			 result <- t
       			 Logic0 <- ((t>=(x-a)) & (t<=(x+a)))  # support Sx={x-a,...,x+a} support de la distribution
			 Logic1 <- ((t<(x-a))|(t>(x+a)))  
       			 tval <- result[Logic0]
			result[Logic1]=0
	     		result[Logic0]<-  ((a+1)^h - (abs(tval-x))^h)/D # Discrete Triangular 				
			return(result)
	}

}



diracDU<-
function(x,t,h,c)
# INPUTS:
#   "c" 	: the number of categories in the Aitchison and Aitken kernel.
#   "x" 	: the target.
#   "t" 	: the single value where the function is computed.
#   "h"	 	: the bandwidth.It is in (0,1] for the binomial kernel.

# OUTPUT:
# Returns the discrete associated kernel estimation at t.
{	
			 result<-t
       			 Logic1 <- (t==x) 
			 Logic0 <- (t!=x)
       			 result[Logic1]<-(1-h)
        		 result[Logic0]<- (h/(c-1)) 
							
        		return(result)
			
}


      
	if(ker=="bino"){	
			result <- t
       			 Logic0 <- (t <= x+1) # support Sx={0,1,...,x+1}
			 Logic1 <- (x+1 < t)
       			 tval <- result[Logic0]
			result[Logic1]=0
        		result[Logic0]<- dbinom(tval,x+1,(x+h)/(x+1))  # The Binomial kernel 
											
        		return(result)
			}
	   

	else	if(ker=="triang"){
		 result <- dtrg(x,t,h,a)  # The discrete Triangular kernel
		return(result)

		}


	  else   if(ker=="dirDU")
		{	
		  result <- diracDU(x,t,h,c)  # The Dirac Discrete Uniform kernel
		   return(result)
		}
    


 }
