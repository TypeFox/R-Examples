gpr = function(logtheta, covfunc.gpr, x, y, xstar = NULL,  partial.derivatives = FALSE)
{
	#toprint = TRUE		
	#if(toprint)print("in gprt to compare with matlab \n")	
	#if(toprint)print(paste(logtheta[1],covfunc.gpr,x[1,1],y[1],xstar[1,1],sep="; "))
	n= dim(x)[1]
	D = dim(x)[2]
	toprint = FALSE		
	tempF1=strsplit(covfunc.gpr, ",")
	covfunc1 = tempF1[[1]][1]
	covfunc23 = paste(tempF1[[1]][2] , tempF1[[1]][3], sep=",")	
	if (length(tempF1[[1]])==3)
	{
		mo= (eval(call( covfunc1, covfunc23)))
		if (mo[1]  != dim(logtheta)[1])  
		{
			print("Error: Number of parameters do not agree with covariance function")
			return (-1)	
		}
		K = eval(call(covfunc1, covfunc23, logtheta ,  x  ))[[1]]
	}else{
			K = eval(call(covfunc.gpr, logtheta, x)) [[1]]    
	}
	L = t(chol(K))	
	alpha = solve(t(L),solve(L,y))    

	if (is.null(xstar))  #nargin ==4
	{
	    	out1 = 0.5*t(y)%*%alpha + sum(log(diag(L))) + 0.5*n*log(2*pi)
	    	out2=0	
		     if( partial.derivatives == TRUE)
		     {         
		    	out2 = rep(0,length(logtheta))       # set the size of the derivative vector		   	 
		   	 W =    solve( t(L) ,  (solve(L, diag(n)) )  )  -  ( alpha %*% t(alpha) )  # precompute for convenience
		    	for (i in 1:length(out2))
		    	{
		    		if (length(tempF1[[1]])==1)
				{
					out2[i] = sum(W*eval(call(covfunc.gpr, logtheta, x, i ,  FALSE)) [[1]] )/2      
		      		}
			      	if (length(tempF1[[1]])==3)
				{
					out2[i] = sum(W*eval(call(covfunc1, covfunc23, logtheta ,  x, i  ,  FALSE))[[1]])/2  

				}
			}
		}# end partial derivative	
	
	}else if(!is.null(xstar))
	{             
  				K = eval(call(covfunc1, covfunc23, logtheta ,  x,  xstar ,TRUE ))
				Kss= K[[1]] 
				Kstar = K[[2]]	  				
				out1 = t(Kstar) %*% alpha                                      # predicted means
				out2 = 0
				if(partial.derivatives == TRUE)
				{
    					if( !is.null(dim(x))){  #can not say if it is null, becasue if it is an array it gives an error
		    					v = solve(L,Kstar)#dim v: 27,51
		    					out2 = Kss - t(colSums(v * v))	    				
    					}
    				}# end partial derivative	
	}
	mout1 =  out1
	mout2 = t(out2)
	result = list(mout1, mout2)	
	return(result)
}
