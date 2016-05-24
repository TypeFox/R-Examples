sq_dist = function(A , B = NULL)
{	
	mym= 0
	resut = 0	
	if( is.null(B)) {
	B=A
	}
	toprint=FALSE
	if(dim(A)[1] != dim( B)[1] ){
		print("Error: column lengths must agree in sq_dist")
	return 
	}
	if(length(A)==1  && length(B)==1){
	A = as.vector(A)
	B =	as.vector(B)
	}
	# result= ((as.matrix(dist(cbind (t(A[]),t(B[])))))^2)/2
	n=dim(A)[2]
	m=dim(B)[2]	 	 
	C= array(0, dim=c(n,m))
		
	if(m==1)	{    #special case, R automatically turns a column matrix  into a one row matrix
		for( d in 1:dim(A)[1]){
			C= C+t((B[rep(d,n),]    -  ( t(A[rep(d,m),]) ))^2 )
		}
	}else{
		for( d in 1:dim(A)[1]){
    			#C = C + (repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2
    			C = C +   (B[rep(d,n),]    -  ( t(A[rep(d,m),]) ))^2
		}
	}
	return (C)
}

covSEiso =function(loghyper = NULL , x = NULL , z =  NULL , testset.covariances= FALSE)
{
	A=B=0
	toprint=FALSE

	if (is.null(loghyper) ) 
	 {
	 	 return(2) 
	}              # report number of parameters
	n = dim(x)[1]
	D = dim(x)[2]
	ell = exp(loghyper[1])                           # characteristic length scale
	sf2 = exp(2*loghyper[2])                                     # signal variance
	A=B= array()
	if (is.null(z)){  
		# as.matrix(dist(cbind(x,y)))
		A = sf2*exp(-sq_dist(t(x)/ell)/2)
		B=0
	}else if (testset.covariances== TRUE )                            # compute test set covariances
    	{
    		if (is.null(dim(z)  ))
    		{
    			matlab.dim = length(z)
    		}else{
    			matlab.dim = dim(z)[1]
    		}
    		A=sf2*rep(1, matlab.dim)
  		B = sf2*exp(-sq_dist(t(x)/ell,t(z)/ell)/2)	
  	}else if (testset.covariances== FALSE ) 
  	{                                                #compute derivative matrix
  		if (length(z)== 1 && z == 1)#dim returns null if there is one number, different from matlab 
  			{       
  			#->                               # first parameter
  			A = sf2*exp(-sq_dist(t(x)/ell)/2)*sq_dist(t(x)/ell) # A = sf2*exp(-sq_dist(x'/ell)/2).*sq_dist(x'/ell);  
	
  			}else{                                       #second parameter
    				A = 2*sf2*exp(-sq_dist(t(x)/ell)/2)
    			}
    	B=  0
 	} 
 	result  = list(A,B)
 	return (result)
}