Conakereport <-
function(Vec,ker,h=NULL,a=0,b=1){
###########################################################################################################
# INPUTS:
#   "Vec" 		: Sample of data 
#   "h" 		: Bandwidth.
#   "ker" 		: The kernel function:  "dirDU" DiracDU,"bino" Binomial,"triang" discrete Triangular. 		   
#   "a" 		: The lef.
#   "b" 		: The rigth bound.
# OUTPUT:Returns a list containing:
#    "h_n" 		: The bandwidth obtained by cross validation.
#    "C_n" 		: The normalizing constant.
###########################################################################################################

 

if(missing(h)){
  		h1=cvbw(Vec,NULL,ker,a,b)
		h=h1$hcv
  		
		 
	bilan=dke(Vec,ker,h,NULL,a,b)
		message('f_n is the estimated p.d.f. obtained')

	if (ker=="GA"){
 		 	message('using gamma kernel and h_n by cross validation technique.') 
	return(list(hn_cv=h,C_n=bilan$C_n))
	}

	else if (ker=="BE") {
			message('using extended beta kernel and h_n by cross validation technique.')
	return(list(hn_cv=h,C_n=bilan$C_n))
 	}

	else if (ker=="LN"){
			message('using lognormal kernel and h_n by cross validation technique.')
                  return(list(hn_cv=h,C_n=bilan$C_n))
		}

	else if (ker=="RIG") {
		message(paste('using reciprocal inverse Gaussian and h_n by cross validation technique.'))
                return(list(hn_cv=h,C_n=bilan$C_n))
		}
    
    
 	}

else 
  {
	bilan=dke(Vec,ker,h,NULL,a,b)
	message('f_n is the estimated p.d.f. obtained')
	if (ker=="GA"){
	message('with a gamma kernel and a given bandwidth h=',h)
        return(list( h=h,C_n=bilan$C_n))
	}
 	
	else if (ker=="BE"){
	message('with an extended beta kernel and a given bandwidth h=',h)    
        return(list( h=h,C_n=bilan$C_n))
	}
 	
	else if (ker=="LN"){
	message('with a lognormal kernel and a given bandwidth h=',h)
       return(list(h=h,C_n=bilan$C_n))
	}
	else if (ker=="RIG"){
	message('with a reciprocal inverse Gaussian kernel and a given bandwidth h=',h)
       return(list(h=h,C_n=bilan$C_n))
	}

 
  }
}
