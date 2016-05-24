	#########################################################################
	#			    Fast ECP (eFast)    			#
	#########################################################################


#Return a list with the following information:
#	number - the number of change points
#	estimates - the change point locations
#	gofM - the gof vector
#	cpLoc - list of all optimal change point locations for differing 
#		numbers of change points
#	time - total running time
#Returned change points indicate the end of a homogeneous segment


#Function arguments are as follows:
#	Z - The time series in which change points are to be found
#	K - The maximum number of possible change points
#	delta - Window size to use when calculating the O(n^2) portion of the U-statistic
#	alpha - Scaling parameter used in distance calculations
#	eps - Epsilon probability used for pruning
#	verbose - Should status updates be printed


e.cp3o = function(Z, K=1, delta=29, alpha=1, eps=0.01, verbose=FALSE){
	#Argument checking
	if(!is.matrix(Z))
		stop("Z must be an n x d matrix.")
	if(alpha <= 0 || alpha > 2)
		stop("alpha must be in the interval (0,2].")
	if(delta < 2)
		stop("delta must be a positive integer greater than 1.")
	if(K < 1 || K > floor(nrow(Z)/(delta+1)))
		stop("K is not in an acceptable range.")
	if(eps <= 0 || eps >= 1)
		stop("eps must be in the interval (0,1).")
	#Force K and delta to be integers
	delta = as.integer(delta)
	K = as.integer(K)

	#Call C++ code that implements the method and store result in res
	#Also keep track of time
	t1 = proc.time()
	res = eFastC(Z, K, delta, alpha, eps, verbose)
	t2 = proc.time()
	res$time = as.numeric((t2-t1)[3])
	#Correct for the fact that C++ is zero based
	#res$estimates = res$estimates + 1
	#res$cpLoc = res$cpLoc + 1

	return(res)
}
