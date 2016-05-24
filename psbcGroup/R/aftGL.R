
aftGL <- function(survData,
					grpInx,
                    hyperParams,
                    startValues,
                    numReps,
                    thin,
                    burninPerc = 0.5
                    )
{	
	###
	n	<- dim(survData)[1]
	p	<- dim(survData)[2] - 2
	
	K <- length(unique(grpInx))

	
	###
	
	nStore <- numReps/thin * (1 - burninPerc)

	mcmc <- .C("aftGLmcmc",
						survData 		= as.double(as.matrix(survData)),
						grpInx			= as.double(grpInx),
						n				= as.integer(n),
						p				= as.integer(p),
						K				= as.integer(K),
                        hyperParams 	= as.double(hyperParams),
                        startValues 	= as.double(startValues),    
						burninPerc      = as.double(burninPerc),
						numReps			= as.integer(numReps),
						thin			= as.integer(thin),
                        samples_alpha   = as.double(rep(0, nStore*1)),
                        samples_beta 	= as.double(rep(0, nStore*p)),
                        samples_sigSq   = as.double(rep(0, nStore*1)),
                        samples_tauSq   = as.double(rep(0, nStore*p)),
                        samples_lambdaSq= as.double(rep(0, nStore*1)),
                        samples_w       = as.double(rep(0, nStore*n))
                        )
    
    alpha.p 	<- matrix(mcmc$samples_alpha, nrow = nStore, byrow = TRUE)
    
	if(p > 0){
		beta.p 		<- matrix(mcmc$samples_beta, nrow = nStore, byrow = TRUE)
	}
	if(p == 0){
		beta.p 		<- NULL
	}
    
    sigSq.p 	<- matrix(mcmc$samples_sigSq, nrow = nStore, byrow = TRUE)
    
	if(p > 0){
		tauSq.p 	<- matrix(mcmc$samples_tauSq, nrow = nStore, byrow = TRUE)
	}
	if(p == 0){
        tauSq.p 	<- NULL
	}
    lambdaSq.p 	<- matrix(mcmc$samples_lambdaSq, nrow = nStore, byrow = TRUE)
    w.p         <- matrix(mcmc$samples_w, nrow = nStore, byrow = TRUE)

    cenInx = which(survData[,2] == 0)
	
	ret <- list(alpha.p = alpha.p, beta.p = beta.p, sigSq.p = sigSq.p, tauSq.p = tauSq.p, lambdaSq.p = lambdaSq.p, w.p = w.p, cenInx = cenInx, data = survData, grpInx = grpInx)

	class(ret) <- "aftGL"
	return(ret)
	
}










