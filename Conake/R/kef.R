kef <-
function(x,t,h,ker,a=0,b=1){
###########################################################################################################
# INPUTS:
#   "x" 	: the target.
#   "t" 	: the single or the grid value where the function is computed.
#   "h" 	: the bandwidth parameter.
#   "ker" 	: the kernel:  "BA" beta extended,"GA" gamma,"LN" lognormal,"RIG" reciprocal inverse Gaussian.	 
#   "a" 	: the left bound of the support of the distribution for extended beta kernel. Default value is 0.
#   "b" 	: the right bound of the support of the distribution for extended beta kernel. Default value is 1.
# OUTPUT:
# Returns the discrete associated kernel value at t.
###########################################################################################################
       
	 if(ker=="BE"){	
			 
			 result <- t
       			 Logic0 <- ((a<=t)&(t<= b)) # support 
			 Logic1 <- ((t<a)|(b<t))
       			 tval <- result[Logic0]
			 result[Logic1]=0
        		
			 result[Logic0]<- ((1/((b-a)^(1+h^(-1))*beta(((x-a)/((b-a)*h))+1,((b-x)/((b-a)*h))+1))))*((tval-a)^((x-a)/((b-a)*h)))*((b-tval)^((b-x)/((b-a)*h))) 
			
        		  return(result)
					
        		}
			
	   
 	  else if(ker=="GA"){
 			result <- t
       			 Logic0 <- (0<=t) # support 
			 Logic1 <- (t<0)
       			 tval <- result[Logic0]
			 result[Logic1]=0
        		 #result[Logic0]<- ((tval^(x/h))/gamma((x/h)+1))*h^(((-x/h)-1))*exp((-tval/h))				
			 result[Logic0]<- dgamma(tval,(x/h)+1,1/h)
        		 return(result)			
			
		 	}

	  else if(ker=="LN"){ 
 			 result <- t
       			 Logic0 <- (0<=t) # support 
			 Logic1 <- (t<0)
       			 tval <- result[Logic0]
			 result[Logic1]=0
        		# result[Logic0]<-  (1/(tval*h*sqrt(2*pi)))*exp((-1/2)*((1/h)*log(tval/x)-h)^2) 				
			 result[Logic0]<- dlnorm(tval,meanlog=log(x)+h^2,sdlog=h)
        		return(result)	     	   
			
		}
	  else if(ker=="RIG"){ 
 			 result <- t
       			 Logic0 <- (0<t) # support 
			 Logic1 <- (t<=0)
       			 tval <- result[Logic0]
			 result[Logic1]<- 0
			 eps<-sqrt(x^2+x*h) # see Libengué (2013)
			 # eps<-1/(x-h)	     # see Scaillet (2013)
			  result[Logic0]<- (1/sqrt(2*pi*h*tval))*exp((-eps/(2*h))*((tval/eps) -2+(eps/tval)))	
        		  return(result)     	   
			
		}
	 
   
 }
