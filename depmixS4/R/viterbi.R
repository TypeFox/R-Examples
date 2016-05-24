# 
# Maarten Speekenbrink, 23-3-2008
# 

viterbi <-
function(object,na.allow=TRUE) {
	# returns the most likely state sequence
	nt <- sum(object@ntimes)
	lt <- length(object@ntimes)
	et <- cumsum(object@ntimes)
	bt <- c(1,et[-lt]+1)
		
	ns <- object@nstates
	
	delta <- psi <- matrix(nrow=nt,ncol=ns)
	state <- vector(length=nt)
	
	prior <- object@init
	
	if(max(ntimes(object)>1)) A <- object@trDens
	B <- object@dens
	if(na.allow) B <- replace(B,is.na(B),1)
	B <- apply(B,c(1,3),prod)
	
	for(case in 1:lt) {
		# initialization
		delta[bt[case],] <- prior[case,]*B[bt[case],]
		delta[bt[case],] <- delta[bt[case],]/(sum(delta[bt[case],]))
		psi[bt[case],] <- 0
		# recursion
		if(object@ntimes[case]>1) {
			for(tt in ((bt[case]+1):et[case])) {
				for(j in 1:ns) {
					if(!object@homogeneous) {
						delta[tt,j] <- max(delta[tt-1,]*(A[tt,j,]))*B[tt,j]
						k <- which.max(delta[tt-1,]*A[tt,j,])
					} else {
						delta[tt,j] <- max(delta[tt-1,]*(A[1,j,]))*B[tt,j]
						k <- which.max(delta[tt-1,]*A[1,j,])
					}
					if(length(k) == 0) k <- 0 # what's this doing here??? can this ever occur? FIX ME
					psi[tt,j] <- k
				}
				delta[tt,] <- delta[tt,]/(sum(delta[tt,]))

			}
		}
		
		# trace maximum likely state
		state[et[case]] <- which.max(delta[et[case],])
		
		# this doesn't need a for loop does it???? FIX ME
		if(object@ntimes[case]>1) {
			for(i in (et[case]-1):bt[case]) {
				state[i] <- psi[i+1,state[i+1]]
			}
		}
	}
  
  colnames(delta) <- paste("S",1:ns,sep="")
  
	delta <- data.frame(state,delta) 	
	return(delta)
}

