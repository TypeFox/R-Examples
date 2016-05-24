#######################
# Function  ess
#   Inputs: 
#			x, a numeric vector of length N assumed to be samples from a Markov chain
#           ignoreBurnin, logical indictating whether or not the first burninProportion of vector x should be ignored
#           burninProportion, if ignoreBurnin == TRUE, the first  burninProportion*length(x)  samples are removed from x before the ess is calculated
#   Outputs:
#           N_ess, the effective number of indepedent samples in the vector
# example:
# 			ess(runif(100))  ## returns a number close to 100
########################

ess <- function(x, ignoreBurnin = FALSE, burninProportion = 0.1){
	if (ignoreBurnin){
		# discard first part of x, according to burninProportion
		if (burninProportion <= 0 || burninProportion >= 1) stop("burninProportion should be in (0,1)")
		x = x[(round(length(x)*burninProportion)+1):length(x)]
	}

	# length of x
	N  = length(x)
	
	if (N < 5) stop("x is too short to get a meaningful ESS")
	 
	# set maximum lag
	maxlag = min(500,N - 2)
	
	r = rep(0,maxlag+1)
	for (k in 0:maxlag){
		r[k+1] = cor(x[1:(N-k)],x[(1+k):N])
	}
	
	G = r[2:(maxlag+1)] + r[1:maxlag]
	
	# estimate autocorrelation
	tauf = -r[1] + 2*G[1]
	for (M in 1:(maxlag-1)){
      if (G[M+1]< G[M] & G[M+1]>0){ 
         tauf = tauf + 2*G[M+1]
      } else {	
      	break
      }
    }	
    tauf = max(tauf,1)
    
    N_ess = N/tauf
    
    return(N_ess)
}
