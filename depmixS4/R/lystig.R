# 
# Ingmar Visser
# 
# LYSTIG algoritme voor de loglikelihood, 23-3-2008
# 

lystig <- function(init,A,B,ntimes=NULL,homogeneous=TRUE,na.allow=TRUE) {

	# Log likelihood computation according to Lystig & Hughes (2002).  This
	# is very similar to the Forward part of the Forward-Backward algorithm
	# but admits of easy computation of the gradients of parameters and the
	# observed information.  This version of the routine only computes the
	# likelihood though.
	
	# NOTE THE CHANGE IN FROM ROW TO COLUMN SUCH THAT TRANSPOSING A IS NOT NECCESSARY ANYMORE 
	# IN COMPUTING ALPHA AND BETA BUT IS NOW NECCESSARY IN COMPUTING XI
	# A = K*K matrix with transition probabilities, from column to row !!!!!!!
	# change to T*K*K
    
	if(na.allow) B <- replace(B,is.na(B),1)
	B <- apply(B,c(1,3),prod)

	# B = T*K*nresp matrix with elements ab_{tij} = P(y_t_i|s_j)
	# init = K vector with initial probabilities !!!OR!!! K*length(ntimes) matrix with initial probs per case
	
	# Returns: 'sca'le factors, recurrent variables 'phi', loglikelihood
		
	nt <- nrow(B)
	ns <- ncol(init)
		
	if(!is.null(ntimes)) {
		phi <- matrix(ncol=ns,nrow=nt)
		sca <- vector(length=nt)
		
		lt <- length(ntimes)
		et <- cumsum(ntimes)
		bt <- c(1,et[-lt]+1)
				
# 		ll <- 0	
		
		for(case in 1:lt) { # multiple cases
			phi[bt[case],] <- init[case,]*B[bt[case],] # initialize case
			sca[bt[case]] <- 1/sum(phi[bt[case],])
			if(ntimes[case]>1) {
				for(i in (bt[case]+1):et[case]) {
					if(homogeneous) phi[i,] <- (A[1,,]%*%phi[i-1,])*B[i,]
					else phi[i,] <- (A[i-1,,]%*%phi[i-1,])*B[i,]
					phi[i,] <- sca[i-1]*phi[i,]
					sca[i] <- 1/sum(phi[i,])
				}
			}
		}
				
		logLike=-sum(log(sca))
		return(list(phi=phi,sca=sca,logLike=logLike))
		
	} else { # single case (single time series, no ntimes provided)
		phi <- matrix(ncol=ns,nrow=nt)
		sca <- vector(length=nt)
		phi[1,] <- init[1,]*B[1,] # initialize
		sca[1] <- 1/sum(phi[1,])
		for(i in 2:nt) {
			if(homogeneous) phi[i,] <- (A[1,,]%*%phi[i-1,])*B[i,]
			else phi[i,] <- (A[i-1,,]%*%phi[i-1,])*B[i,]			
			phi[i,] <- sca[i-1]*phi[i,]
			sca[i] <- 1/sum(phi[i,])
		}
		
		logLike=-sum(log(sca))
		return(list(phi=phi,sca=sca,logLike=logLike))
	}
}

