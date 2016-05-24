

#########################################################################
# compute log-likelihood for din objects
vcov.loglike.din <- function( weights , skillprobs0 , slip0 , guess0 ,
	     latresp , item.patt.split , resp.ind.list ,
		 return.p.xi.aj=FALSE ){
	########################
	IP <- N <- length(weights)
	L <- length(skillprobs0)
	J <- length(guess0)
	# calculate probabilities
	slipM <- matrix( slip0 , nrow= nrow(latresp) , ncol=ncol(latresp))
	guessM <- matrix( guess0 , nrow= nrow(latresp) , ncol=ncol(latresp))
	pj <- (1 - slipM )*latresp + guessM * ( 1 - latresp )
	pjM <- array( NA , dim=c(J,2,L) )
	pjM[,1,] <- 1 - pj
	pjM[,2,] <- pj
	skillprobsM <- matrix( skillprobs0 , nrow=IP , ncol=L , byrow=TRUE )
	# calculate log-likelihood
	h1 <- matrix( 1 , nrow=IP , ncol=L )
	res.hwt <- calc_posterior.v2(rprobs= pjM , gwt=h1 , resp=item.patt.split , 
									 nitems= J , 
									 resp.ind.list=resp.ind.list , normalization=FALSE , 
									 thetasamp.density= NULL , snodes=0 )
	p.xi.aj <- res.hwt$hwt
	# Log-Likelihood (casewise)
	ll2 <- log( rowSums( p.xi.aj * skillprobsM ) )
	if (return.p.xi.aj){	
		res <- list( "ll" = ll2 , "p.xi.aj"=p.xi.aj )			
			}  else { 
		res <- ll2 
				}
	return(res)
			}
#########################################################################	