######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: covNoise
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
covNoise = function(loghyper= NULL , x = NULL , z = NULL ,  testset.covariances= FALSE){	
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
######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: sq_dist
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
sq_dist = function(A , B = NULL){	
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
######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: covSEiso
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
covSEiso =function(loghyper = NULL , x = NULL , z =  NULL , testset.covariances= FALSE){
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
######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: covSum
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
covSum= function(covfuncsum , logtheta  = NULL,  x  = NULL, z = NULL, testset.covariances= FALSE){
	
	A=B=0
	toprint= FALSE	
	temp1=strsplit(covfuncsum , ",")
	covfuncs1 =temp1[[1]][1]
	covfuncs2 =temp1[[1]][2]
	covarray = as.vector(c(covfuncs1, covfuncs2 ))	
	j=as.array( c(eval(call(covfuncs1) ) , eval(call(covfuncs2 ))) )	#j stores number of parameters: j=[2 1]
	v = NULL
	for (i in 1:length(j)){
		v= cbind( v, array(rep(i,j[i] ), dim= c(1,j[i]))		)
	}
	#v is parameter mapper array. 
	#v= [1,1,2]  ,means for the fist 2 parameter of loghyper must be taken for the fist function
	# and third parameter of loghyper must be used for the second function
	if (is.null(logtheta))
	{                                  # report number of parameters
	  A =   eval(call(covfuncs1) ) + eval(call(covfuncs2 ))
	  B=0
  
	}else{
		n = dim(x)[1]
		D = dim(x)[2]
	#	v = array()            # v vector indicates to which covariance parameters belong


	  	if(is.null(z)){                        #nargin case==3                      # compute covariance matrix
	   		A1 = eval(call(covfuncs1, logtheta[v==1], x)) [[1]]
	   		A2= eval(call(covfuncs2, logtheta[v==2], x))[[1]]
	   		A=A1+A2	
	   		B=0	  
			}else{		#nargin case==4
				if (testset.covariances == TRUE){ #nargout==2
			    		ans1 =    eval(call(covfuncs1, logtheta[v==1], x, z, testset.covariances))             
			  	 	ans2=    eval(call(covfuncs2, logtheta[v==2], x, z, testset.covariances))
			      		A = ans1[[1]] [1]+ ans2[[1]]			      	 
					B = ans1[[2] ]+ans2[[2]]     	 	  
				}else if (testset.covariances == FALSE){
					i=v[z]
					j=sum(v[1:z]==i)
					f = covarray[i]
#print(paste("logtheta[v==i] ,f,z,i,j:",logtheta[v==i] ,f,z,i,j,sep="; "))
					A= eval(call(f, logtheta[v==i], x, j, testset.covariances)) [[1]]
					B=0		
				}
			}
		}	
		result = list(A,B)
		return( result)
}
######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: gpr
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
gpr = function(logtheta, covfunc.gpr, x, y, xstar = NULL,  partial.derivatives = FALSE){
	#toprint = TRUE		
	#if(toprint)print("in gprt to compare with matlab \n")	
	#if(toprint)print(paste(logtheta[1],covfunc.gpr,x[1,1],y[1],xstar[1,1],sep="; "))
	
	n= dim(x)[1]
	D = dim(x)[2]
	
	tempF1 = strsplit(covfunc.gpr, ",")
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
	}
	else
	{
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
			W = solve( t(L) ,  (solve(L, diag(n)) )  )  -  ( alpha %*% t(alpha) )  # precompute for convenience
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
	}
	else if(!is.null(xstar))
	{             
		K = eval(call(covfunc1, covfunc23, logtheta ,  x,  xstar ,TRUE ))
		Kss= K[[1]] 
		Kstar = K[[2]]	  				
		out1 = t(Kstar) %*% alpha                                      # predicted means
		out2 = 0
		if(partial.derivatives == TRUE)
		{
			if( !is.null(dim(x)))
			{  #can not say if it is null, because if it is an array it gives an error
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
######################################################################################################################
#																													 
#                     Software that implements										
# 
#          GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION
#
# Copyright (c) 2005 - 2007 by Carl Edward Rasmussen and Chris Williams
#
# Permission is granted for anyone to copy, use, or modify these programs for
# purposes of research or education, provided this copyright notice is retained,
# and note is made of any changes that have been made.
#
# These programs are distributed without any warranty, express or
# implied. As these programs were written for research purposes only, they
# have not been tested to the degree that would be advisable in any
# important application.  All use of these programs is entirely at the
# user's own risk.
#
# The code and associated documentation are available from
#
# http://www.GaussianProcess.org/gpml/code
#
######################################################################################################################
#
#	Function: minimize
#	This function is an R version of MATLAB implementation (gpml) by Carl Edward Rasmussen and Chris Williams.
#	This R version is derived from "gpr" package (http://cran.r-project.org/web/packages/gpr/)
#
######################################################################################################################
minimize = function (X, f, .length, covfunc, x, y){
	INT = 0.1
	EXT = 3.0                
	MAX = 20                  
	RATIO = 10             
	SIG = 0.1
	RHO = SIG/2	 
	if(is.array( .length)){
		if (max(dim( .length)) == 2)
		{
			 red= .length[2] 
			 .length=  .length[1]
		}
	} else{ 
			red=1
	}
	if ( .length>0)
	{
		 S='Linesearch'
	} else{
		 S='Function evaluation'
	} 
	
	i.m= 0                                            # zero the run length counter
	ls_failed = 0                            # no previous line search has failed
	f_out =eval(call (f, X, covfunc , x, y, NULL , TRUE))
	#print(paste("fout", f_out,sep="|"));
	f0= f_out[1][[1]]
	df0 = t(f_out[2][[1]]) #out put is a colum, R makes it a row when put it in list to output, I conver it back here
	
	
	fX = f0
	i.m = i.m + ( .length<0)                                            # count epochs?!
	
	s = -df0
	d0 = -t(s)%*%s           # initial search direction (steepest) and slope
	x3 = red/(1-d0)                      # initial step is red/(|s|+1)
	
	mainloop = TRUE
	while( i.m < abs( .length) && mainloop)
	{                                      # while not finished
  		i.m = i.m + (.length > 0)                                     # count iterations?!
  		X0 = X
  		F0 = f0
  		dF0 = df0                  # make a copy of current values
		
  		if ( .length>0)
			M = MAX
  		else
			M = min(MAX, - .length- i.m)
			
		
			
  		whilerun = TRUE
		f3_=c(0)
		df3_=c(0)
 		while (whilerun ==TRUE)                             # keep extrapolating as long as necessary
    		{
    			x2 = 0
    			f2 = f0
    			d2 = d0
    			f3 = f0
    			df3 = df0
    			success = FALSE
    			while  (success == FALSE && M > 0)
    			{
    				M = M - 1
    				i.m = i.m + (.length<0)                         #count epochs?!
    				options(show.error.messages = FALSE)
					f_out2 =eval(call (f , X+ x3[1]*s, covfunc , x, y, NULL, TRUE))
					
					f3=  f_out2[1][[ 1 ]][[ 1 ]]
					df3 = t(f_out2[2][[1]])
					
					f3_=rbind(f3_,f3,0);
					df3_=rbind(df3_,df3,0);
					
					if (is.na(f3) || is.infinite(f3) || is.nan(f3)  || any(is.nan(df3) || is.na(df3) || is.infinite(df3))  )  
					{
						  cat( " ")
						  x3 = (x2+x3)/2 
					}
					else
					{
						success = TRUE
					}
					options(show.error.messages = TRUE)		  
				}
				
				if (f3 < F0)
				{
					X0 = X+x3[1]*s
					F0 = f3
					dF0 = df3
				}         # keep best values
    			d3 = t(df3)%*%s                                       # new slope
    			if( d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  )# are we done extrapolating?
     			{
     				 whilerun = FALSE
     				 break	 
    			}	
      			x1 = x2
      			f1 = f2
      			d1 = d2                       # move point 2 to point 1
    			x2 = x3 
    			f2 = f3
    			d2 = d3                        # move point 3 to point 2
    			A = 6*(f1-f2)+3*(d2+d1)*(x2-x1)                 # make cubic extrapolation
    			B = 3*(f2-f1)-(2*d1+d2)*(x2-x1)
    			x3 = x1-d1*(x2-x1)^2/(B+sqrt(abs(B*B-A*d1*(x2-x1)))) # num. error possible, ok!		
    			if ( (B*B-A*d1*(x2-x1)   < 0)[1]  || is.nan(x3) || is.infinite(x3) || x3 < 0) # num prob | wrong sign?
      			{	
      				x3 = x2*EXT                                 # extrapolate maximum amount
    			}else if( x3 > x2*EXT)                  # new point beyond extrapolation limit?
    			{
      				x3 = x2*EXT                                 # extrapolate maximum amount
      			}else if (x3 < x2+INT*(x2-x1) )         # new point too close to previous point?
      			{	
      				x3 = x2+INT*(x2-x1)
   			}
   			#x3= round(x3*10000)/10000   #this is just to make same answers with matlab code
      		}
      		#print(f3_)
			#print(df3_)
			#ANSWER <- readline("continue? ")
			#if (substr(ANSWER, 1, 1) == "n")
			#	stop();
			
			while ((abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0 ) # keep interpolating
      		{
    			if( d3 > 0 || f3 > f0+x3*RHO*d0)                         # choose subinterval
    			{
      				x4 = x3
      				f4 = f3
      				d4 = d3                      #move point 3 to point 4
      			}else{
      				x2 = x3
      				f2 = f3
      				d2 = d3                      # move point 3 to point 2
      			}
    			if (f4 > f0 ){          
      				x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2))  # quadratic interpolation
    			}else{
      				A = 6*(f2-f4)/(x4-x2)+3*(d4+d2)                    # cubic interpolation
      				B = 3*(f4-f2)-(2*d2+d4)*(x4-x2)
      				x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A        # num. error possible, ok!
    			}
      			
      			 if (is.nan(x3) || is.infinite(x3) )
      			 {
      				x3 = (x2+x4)/2               # if we had a numerical problem then bisect
    			 }
    			x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2))  # don't accept too close
    			f_out3 =eval(call (f , X+ x3*s, covfunc , x, y, NULL, TRUE))  
			f3 = f_out3[1] [[1]][[1]]
			df3 = t(f_out3[2][[1]])
    			if (f3 < F0)
    			{
    				x3=x3[[1]] # arrays of one in matlab converst to number, in R it is not
    				 X0 = X+x3*s
    				F0 = f3
    				dF0 = df3
    			}# keep best values
    			M = M - 1
    			i.m = i.m + (.length<0)                            # count epochs?!
    			d3 =  t(df3)%*%s                                                    # new slope
            }#end while                                        #end interpolation
 	   if (abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0 )        # if line search succeeded
 	   {
 	   	x3=x3[[1]]
    		X = X+x3*s
    		f0 = f3
    		fX=	t(cbind (t(fX), f0))
    		 s=(  (   (t(df3)%*%df3- t(df0)%*%(df3))[1]  )/ (  (t(df0)%*%(df0))[[1]]    ) *s  ) - df3
    		df0 = df3                                               # swap derivatives
    		d3 = d0
    		d0 = t(df0)%*%s
    		if (d0 > 0 ) 
    		{                                    # new slope must be negative
      			s = -df0
      			d0 = -t(s)%*%s                  # otherwise use steepest direction
    		}
    		x3 = x3 * min(RATIO, d3/(d0-  (  2^(-1022) )  ))          # slope ratio but max RATIO
   		 ls_failed = 0                           # this line search did not fail
 	  } else{
    		X = X0
    		f0 = F0
    		df0 = dF0                     # restore best point so far
   		 if (ls_failed || i.m > abs(.length)   )    # line search failed twice in a row
    		{
      			mainloop = 0 #break; 	                            #or we ran out of time, so we give up
      			break
    		}
    		s = -df0
     		d0 = -t(s)%*%s                                        #try steepest
    		x3 = 1/(1-d0)                    
    		ls_failed = 1                                    # this line search failed
  	  }#end if
	}#end first while
	
  return (list(X))  #i.m is number of epochs
}
########################################################
#													   #
#	Function: normalize								   #	
#	Author: Mohsen Ahmadi							   #
#	Email: mohsen_ahmadi989@yahoo.com  				   #
#	Date: Feb 2014								   	   #
#													   #
########################################################
normalize <- function(x){

	mu = colMeans(x);
	x_norm = sweep(x, length(dim(x)) , mu, FUN="-")
	
	if (class(x) == "data.frame")
	{
		sigma = sapply(x_norm, sd)
	}
	else
	{
		sigma = sd(x_norm)
	}
	
	sigma[sigma==0]=1	  	
	norm = sweep(x_norm, length(dim(x)), sigma, FUN="/" )
	
	return (list(norm, mu, sigma))
}
########################################################
#													   #
#	Function: one_nn_tc								   #	
#	Author: Mohsen Ahmadi							   #
#	Email: mohsen_ahmadi989@yahoo.com  				   #
#	Date: Feb 2014								   	   #
#													   #
########################################################
one_nn_tc = function(X1,X2){
	index=0
	tc=-1
    for (k in 1:dim(X2)[1])
	{
        K=tanimoto_coef(X1,X2[k,])
        #CheckTc(K, X1, X2[k,])
		if ( K>tc)
		{
			tc=K
			index=k
		}
    }
	
    return( index )
}
########################################################
#													   #
#	Function: tanimoto_coef							   #	
#	Author: Mohsen Ahmadi							   #
#	Email: mohsen_ahmadi989@yahoo.com  				   #
#	Date: Feb 2014								   	   #
#													   #
########################################################
tanimoto_coef = function(X, Y){
	
	X = as.matrix(X)#R automatically converts matrixes to data frame, we have to convert them back
	Y = as.matrix(Y)
	res = ( t(X) %*% Y ) / (( t(X) %*% X ) + ( t(Y) %*% Y ) - ( t(X) %*% Y ))
	return(res)
}
########################################################
#													   #
#	Function: SimuChemPC							   #	
#	Author: Mohsen Ahmadi							   #
#	Email: mohsen_ahmadi989@yahoo.com  				   #
#	Date: Feb 2014								   	   #
#													   #
########################################################
SimuChemPC = function(dataX, dataY, method = "RA", experiment = 1){
	options(warn=-1);
	if (is.null(dataX) || is.null(dataY))
	{
	 	  print("Error : Input parameteres can not be null.")
	 	  return
	}
	
	simulation=0
	if (method == "RA") simulation=1
	if (method == "EI") simulation=2
	if (method == "NN") simulation=3
	if (method == "GP") simulation=4
	if(simulation == 0)
	{
		print("Error : method type should be one of: EI, GP, NN or RA")
		return
	}	
	
	seedpath = "seeds_for_random_generatorMatlab.txt"
	seedFile = system.file("extdata", seedpath , package="SimuChemPC")
	a <- read.table(seedFile)   # read data from file with tab delimiter
	Seed = a$V1 #convert the frame into list of integer values	
	seed_i=1;
	len=dim(dataX)[1]
	len_half = floor(len/2);
	potency_real=matrix(data=-1000,nrow=len-len_half,ncol=experiment);
	for(loop in 1:experiment)
	{
		START=1;
        END=len;
        xTest=dataX;
        yTest=log10(dataY);
		
		xTrain=c();
		yTrain=c();
		for(i in 1:len_half)
		{
			ran = round(START + (END-START) * Seed[seed_i])
			seed_i= seed_i+1
			
			xTrain = rbind(xTrain, xTest[ran,]);
            xTest=xTest[-ran,];
            
            yTrain = rbind(yTrain, yTest[ran]);
            yTest=yTest[-ran];
			
			END=END-1
		}
		yTest=as.matrix(yTest);
		#===========================================================
		xNorm = normalize(xTrain);
		yNorm = normalize(yTrain);
		
		x = xNorm[[1]];
		y = yNorm[[1]];

		testX_norm = sweep(xTest, length(dim(xTest)) , xNorm[[2]], FUN="-")
		xstar = sweep(testX_norm, length(dim(xTest)), xNorm[[3]], FUN="/" )
		
		testY_norm = sweep(yTest, length(dim(yTest)) , yNorm[[2]], FUN="-")
		ystar = sweep(testY_norm, length(dim(yTest)), yNorm[[3]], FUN="/" )
		#===========================================================
		p_values = c();
		for(j in 1:dim(x)[2])
		{
			p = cor.test(as.matrix(x[, j]), y, method='spearman')$p.value
			p_values = rbind(p_values, p);
		}
		
		adj_pVals = p.adjust(p_values, method = "BH");
		hVals = as.integer(adj_pVals <= 0.05);
		
		if (sum(hVals) > 0)
		{
			attr = grep(1, hVals);
		}
		else
		{
			min_p = min(adj_pVals);
			attr = grep(TRUE, adj_pVals==min_p);
		}
		x = as.matrix(x[, c(attr)]);
		xstar = as.matrix(xstar[, c(attr)]);
		#===========================================================
		
		counter=1;
		print(paste("experiment => " , loop, sep=""))
		while (TRUE)
		{
			d1=dim(xstar)[1]
			if (simulation==1) # RA
			{  
				index = round(runif(1, 1, d1-1));
			} 
			else #  EI & NN & GP
			{ 
				if (simulation==3) 
				{
					row = which.max(y);
					index = one_nn_tc(x[row,], xstar)
				}
				else #EI & GP
				{ 
					covfunc="covSum,covSEiso,covNoise" 
					loghyper= array(c(-1,-1,-1), dim=c(3,1))
					loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y)
					loghyper = loghyper[[1]];
					f.out= gpr(loghyper, covfunc, x, y, xstar , TRUE)
					mu = f.out[1] [[1]]
					S2 = f.out[2] [[1]]
					if (simulation==4 )
					{
						index = which.max(mu);
					} 
					else 
					{ # EI
						Qmax = max(y)
						ei = rep(1, d1)
						for (j in 1:d1)
						{
							m = mu[j]
							var = S2[j]
							
							index_=one_nn_tc(xstar[j,],x)
							s=sqrt(max( var, (y[ index_] - m ) ^ 2 ))
							
							u = m - Qmax
							ei[j] = (u * pnorm(u/s,0,1)) + (s * dnorm(u/s,0,1))
						}
						index = which.max(ei);
					}
				}
			}
			
			potency_real[counter,loop]=ystar[index];    
			counter = counter + 1;
			
			x = rbind(x, xstar[index,])  
			y = rbind(y, ystar[index])
			
			# remove selected compound
            xstar=xstar[-index,]
			xstar=as.matrix(xstar);
			ystar=ystar[-index]
			
			if (counter==len_half)
			{
				print("=============")
				break;
			}
		} # end while
	} # end for
	#options("scipen"=100, "digits"=6)
	return(potency_real)
}