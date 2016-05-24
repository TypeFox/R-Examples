lowfreqdrift <-
function(dim, freq=128, nscan, TR, template, verbose=TRUE){

	spm_drift <- function(N, K) {
    	  n <- 0:(N-1)
    	  C <- matrix(0, nrow=N, ncol=K)

     	  C[,1] <- 1/sqrt(N)
    	  for(k in 2:K) {
       	  C[,k] = sqrt(2/N)*10*cos(pi*(2*n+1)*(k-1)/(2*N))
    	  }
	  return(C)
	}

	n <- floor( 2*(nscan*TR)/freq + 1 )
	if(n<3){
		stop("Number of basis functions is too low. Lower frequency or longer scanning time should be provided.")
	}
	drift.base <- spm_drift(nscan, n)[,-1]  
        drift.image <- array(rep(1, prod(dim)), dim=c(dim))
  	drift.array <- drift.image %o% rowSums(drift.base)

        if(!missing(template)){
 		if(length(dim(template))>3){
			stop("Template should be a 2D or 3D array.")
		}
                template.time <- array(rep(template,nscan), dim=c(dim,nscan))
                ix <- which(template.time!=0)
                drift.array[-ix] <- 0
        }

	return(drift.array)
}

