covSum= function(covfuncsum , logtheta  = NULL,  x  = NULL, z = NULL, testset.covariances= FALSE)
{
	
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