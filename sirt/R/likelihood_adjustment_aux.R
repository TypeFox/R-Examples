

#######################################################
# compute adjusted likelihood
likelihood_adjustment_compute <- function( likelihood, theta , 
			thetaM , adjfac , tuningfac=1){
	#*******************
	TP <- length(theta)	
	# compute weighted mean and weighted SD
	res0 <- likelihood_moments( likelihood , theta=theta )	
	M1 <- res0$M
	SD1 <- res0$SD
	w1 <- rowSums(likelihood)	
	# compute adjusted likelihood
	like2 <- 0*likelihood
	for (tt in 1:TP){
		like2[,tt] <- stats::dnorm( theta[tt] , mean=M1 , sd = SD1*adjfac*tuningfac )
					}
    like2 <- like2 / rowSums(like2) * w1
    return(like2)
			}
#############################################################			
# EAP reliability
like_adj_EAP_reliability <- function( M , SD ){
	var( M ) / ( var(M) + mean( SD^2 ) )
				}
#########################################################				
# likelihood adjustment with a tuning parameter
likelihood_adjustment_tuning <- function( likelihood , theta , thetaM , adjfac , 
		tuningfac , probs , maxiter=100 , conv=5*1E-5 , normal_approx=TRUE ){
		
	like2 <- likelihood_adjustment_compute( likelihood, theta , thetaM , adjfac ,
				tuningfac=tuningfac)
	change <- 1
	iter <- 0
	
	if (normal_approx){
		probs <- trait_normal_approx( probs , theta )
					 }
	
	
	while( ( iter < maxiter ) & ( change > conv ) ){
	
		probsM <- matrix( probs , nrow(likelihood) , ncol(likelihood) , byrow=TRUE )
		post <- like2 * probsM
		post <- post / rowSums(post)
		probs_new <- colMeans( post )
		if (normal_approx){
			probs_new <- trait_normal_approx( probs_new , theta )
						 }
		
		change <- max( abs( probs_new - probs ))
		probs <- probs_new		
		iter <- iter + 1
					}	
				
	res0 <- likelihood_moments( likelihood = like2 * probsM , theta=theta  )	
	EAP.rel <- like_adj_EAP_reliability( res0$M , res0$SD )
	res <- list( "likelihood" = like2 , "EAP.rel" = EAP.rel )
	return(res)
			}
##########################################################
# normal approximation of trait distribution
trait_normal_approx <- function( probs , theta ){
	M <- sum( theta*probs )
	SD <- sqrt( sum( theta^2*probs ) - M^2 )
	probs1 <- stats::dnorm( theta , mean=M , sd=SD)
	probs1 <- probs1/sum(probs1)			
	return(probs1)
			}
############################################################			