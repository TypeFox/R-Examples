kef <-
function(x,t,h,type_data=c("discrete","continuous"),ker=c("bino","triang",
                                  "dirDU","BE","GA","LN", 
                                  "RIG"),a0=0,a1=1,a=1,c=2){
###########################################################################################################
# INPUTS:
#   "x" 	: the target.
#   "t" 	: the single or the grid value where the function is computed.
#   "h" 	: the bandwidth parameter.
#   "ker" 	: the kernel:  "bino" binomial,"triang" discrete triangular,"dirDU" Dirac discrete uniform.
#			       "BE"  extended beta,"GA" gamma,"LN" lognormal,"RIG" reciprocal inverse Gaussian.     	 
#   "a0" 	: the left bound of the support of the distribution for extended beta kernel. Default value is 0.
#   "a1" 	: the right bound of the support of the distribution for extended beta kernel. Default value is 1.
#   "a" 	: The arm is used only for the discrete triangular distribution.
#   "c" 	: The number of categories in DiracDU kernel and is used only for DiracDU
# OUTPUT:
# Returns the discrete associated kernel value at t.
###########################################################################################################
    if (missing(type_data))  stop("argument 'type_data' is omitted")
    if ((type_data=="discrete") & (ker=="GA"||ker=="LN"||ker=="BE" ||ker=="RIG")) 
       stop(" Not appropriate kernel for type_data")
   if ((type_data=="continuous") & (ker=="bino"||ker=="triang"||ker=="dirDU")) 
      stop(" Not appropriate kernel for 'type_data'")
   if ((type_data=="discrete") & missing(ker)) ker<-"bino"
   if ((type_data=="continuous") & missing(ker)) ker<-"GA"
      
dtrg<-function(x,t,h,a){

if (a==0)
	{
			 result <- t
       			 Logic1 <- (t==x) 
			 Logic0 <- (t!=x) 
       			 result[Logic1]=1
        		 result[Logic0]=0					
        		
			}
  
else

      {	
		
			 u=0:a;
 			 u=sum(u^h)			 
			 D=(2*a+1)*(a+1)^h -2*u                 
			 result <- t
       			 Logic0 <- ((t>=(x-a)) & (t<=(x+a)))  # support Sx={x-a,...,x+a} support of the distribution
			 Logic1 <- ((t<(x-a))|(t>(x+a)))  
       			 tval <- result[Logic0]
			result[Logic1]=0
	     		result[Logic0]<-  ((a+1)^h - (abs(tval-x))^h)/D # Discrete Triangular 				
			
	}
return(result)
}



diracDU<-
function(x,t,h,c)
{	
			 result<-t
       			 Logic1 <- (t==x) 
			 Logic0 <- (t!=x)
       			 result[Logic1]<-(1-h)
        		 result[Logic0]<- (h/(c-1)) 
							
        	
			
}


      
	if(ker=="bino"){	
			result <- t
       			 Logic0 <- (t <= x+1) # support Sx={0,1,...,x+1}
			 Logic1 <- (x+1 < t)
       			 tval <- result[Logic0]
			result[Logic1]=0
        		result[Logic0]<- dbinom(tval,x+1,(x+h)/(x+1))  # The Binomial kernel 
											
        	
			}
	   

	else	if(ker=="triang"){
		 result <- dtrg(x,t,h,a)  # The discrete Triangular kernel
	

		}


	  else   if(ker=="dirDU")
		{	
		  result <- diracDU(x,t,h,c)  # The Dirac Discrete Uniform kernel
		
		}
    


	 if(ker=="BE"){	
			 
			 result <- t
       			 Logic0 <- ((a0<=t)&(t<= a1)) # support 
			 Logic1 <- ((t<a0)|(a1<t))
       			 tval <- result[Logic0]
			 result[Logic1]=0
        		
			 result[Logic0]<- ((1/((a1-a0)^(1+h^(-1))*beta(((x-a0)/((a1-a0)*h))+1,((a1-x)/((a1-a0)*h))+1))))*((tval-a0)^((x-a0)/((a1-a0)*h)))*((a1-tval)^((a1-x)/((a1-a0)*h))) 
			
        		
					
        		}
			
	   
 	  else if(ker=="GA"){
 			result <- t
       			 Logic0 <- (0<=t) # support 
			 Logic1 <- (t<0)
       			 tval <- result[Logic0]
			 result[Logic1]=0		
			 result[Logic0]<- dgamma(tval,(x/h)+1,1/h)
        					
			
		 	}

	  else if(ker=="LN"){ 
 			 result <- t
       			 Logic0 <- (0<=t) # support 
			 Logic1 <- (t<0)
       			 tval <- result[Logic0]
			 result[Logic1]=0
        		 # result[Logic0]<-  (1/(tval*h*sqrt(2*pi)))*exp((-1/2)*((1/h)*log(tval/x)-h)^2) 				
			 result[Logic0]<- dlnorm(tval,meanlog=log(x)+h^2,sdlog=h)
        	     	   
			
		}
	  else if(ker=="RIG"){ 
 			 result <- t
       			 Logic0 <- (0<t) # support 
			 Logic1 <- (t<=0)
       			 tval <- result[Logic0]
			 result[Logic1]<- 0
			 eps<-sqrt(x^2+x*h) # see Libengué (2013)
			 result[Logic0]<- (1/sqrt(2*pi*h*tval))*exp((-eps/(2*h))*((tval/eps) -2+(eps/tval)))	
        		    	   
			
		}
	 
   return(result)
 }
