covNoise = function(loghyper= NULL , x = NULL , z = NULL ,  testset.covariances= FALSE)
{	
	toprint=FALSE
	if (is.null(loghyper) ) 
	 {
	 	 return(1) 
	}              # report number of parameters

		s2 = exp(2*loghyper)                                     # noise variance

	if (is.null(z))
	{                                     # compute covariance matrix
  		A = s2* diag(dim(x)[1])  #I change the code as here it was giving error by %*% (in source code too)
  		B=0
	}else if (testset.covariances== TRUE  )                            # compute test set covariances
  	{
  		A = s2
  		B = 0                              # zeros cross covariance by independence
	}else if (testset.covariances== FALSE  ) 
	{                                                 # compute derivative matrix
  		A = 2*s2 * diag(dim(x)[1])
  		B= 0	
	}
	 result= list(A, B) 	 
	return (result)
}
