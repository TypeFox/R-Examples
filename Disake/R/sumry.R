sumry <-
function(Vec,type_bw,ker="bino",h=NULL,a=1,c=2){
###########################################################################################################
# INPUTS:
#   "Vec" 		: Sample of data 
#   "h" 		: Bandwidth.
#   "ker" 		: The kernel function:  "dirDU" DiracDU,"bino" Binomial,"triang" discrete Triangular. 		   
#   "a" 		: The arm is used only for the Discrete Triangular kernel. The default value is 1.
#   "c" 		: The number of categories in the Aitchison and Aitken kernel  is used only for DiracaDU.The default value is 2.
# OUTPUT:Returns a list containing:
#    "n" 		: The number of observations.
#    "support_of_f_n" 	: The support of f_n.
#    "C_n" 		: The normalizant constant.
#    "ISE_0" 		: The integrated squared error when using the naive distribution instead of f_n.
#    "f_0" 		: The couples (x,f_0(x)).
#    "f_n" 		: The couples (x,f_n(x)).
#    "h" 		: The bandwidth used to estimate the p.m.f.
###########################################################################################################

 

if(missing(h)){
  	if((ker=="dirDU")&(type_bw=="CV")){
		h1=CVbw(Vec,NULL,ker,a,c)
		h=h1$hcv
		
   		}
	else if((ker=="triang")& (type_bw=="CV")){
		h1=CVbw(Vec,NULL,ker,a)
		h=h1$hcv
        	}

     	else if((ker=="bino")&(type_bw=="CV")){
			h1=CVbw(Vec,NULL,ker)
			h=h1$hcv
         		}
        else if((ker=="bino")&(type_bw=="Bays")){
	 		h=Baysbw(Vec)
  			}
  
		 

	bilan=kpmfe(Vec,h,ker,a,c)
		message('The estimated p.m.f. f_n is the smoothing of the empirical p.m.f. f_0')

	if ((type_bw=="CV")&(ker=="bino")){
 		 	message('using Binomial kernel and h_n by cross validation technique.') 
return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, hn_cv=h))
	}

	else if ((type_bw=="CV")&(ker=="triang")) {
			message(paste('using Discrete Triangular kernel with a=', a ,'and h_n by cross validation technique. ', sep=" "))
return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, hn_cv=h))
 	}

	else if((ker=="bino")&(type_bw=="Bays")){
			message('using Binomial kernel and h_n by local Bayesian procedure.')
                  return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, hn_Bays=h))
		}

	else if ((type_bw=="CV")&(ker=="dirDU")) {
		message(paste('using Dirac Discrete Uniform kernel with c=', c ,' and h_n by cross validation technique.'))
                return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, hn_cv=h))
		}
    
 	}

else 
  {
	bilan=kpmfe(Vec,h,ker,a,c)
	message('The estimated p.m.f. f_n is the smoothing of the empirical p.m.f. f_0')
	if (ker=="bino"){
	message('with a Binomial kernel and a given bandwidth h=',h)
        return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, h=h))
	}
 	
	else if (ker=="triang"){
	message('with a Discrete Triangular kernel with a= ' , a , ' and a given bandwidth h=',h)
        return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, h=h))
	}
 	
	else if (ker=="dirDU"){
	message('with a Dirac Discrete Uniform kernel with c=', c ,'  and a given bandwidth h=',h)
       return(list(n=bilan$n,support_of_f_n=bilan$support,f_0=bilan$f_0,f_n=bilan$f_n,ISE_0=bilan$ISE_0,C_n=bilan$C_n, h=h))
	}


 
  }
}
