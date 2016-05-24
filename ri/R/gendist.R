gendist <-
function(Ys,perms,X=NULL,Ypre=NULL,prob=NULL,HT=FALSE) {

	numiter <- ncol(perms)
	distout <- rep(NA,numiter)
	
	if (is.null(prob)) {
		prob <- genprob(perms)
		warning("Generating probabilities from permutation matrix.")
	}

	# switch to apply?

	for (iter in 1:numiter) {
		Zri <- perms[,iter]
		Yri <- Ys$Y0
		Yri[Zri==1] <- Ys$Y1[Zri==1]

		distout[iter] <- estate(Yri,Zri,X,Ypre,prob,HT)
		}		
		
	return(distout)
	}
