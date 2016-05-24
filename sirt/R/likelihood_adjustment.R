

####################################################################
# likelihood adjustment
likelihood.adjustment <- function( likelihood , theta = NULL , prob.theta=NULL ,
			adjfac = rep(1,nrow(likelihood)) , extreme.item = 5 ,
			target.EAP.rel = NULL , min_tuning=.2 , max_tuning = 3 ,
			maxiter = 100 , conv = .0001 , trait.normal=TRUE ){

	like0 <- likelihood
	eps <- 1E-30
	normal_approx <- trait.normal
	
	if ( is.null(theta) ){
		theta <- attr( likelihood , "theta" )[,1]
						}

	if ( is.null(prob.theta) ){
		prob.theta <- attr( likelihood , "prob.theta" )
						}

						
	attr(likelihood,"prob.theta") <- NULL
	attr(likelihood,"theta") <- NULL
	attr(likelihood,"G") <- NULL
	
	#**********************
	# add extreme item
	N <- nrow(like0)
	TP <- length(theta)
	thetaM <- matrix( theta , nrow=N , ncol=TP  , byrow=TRUE)
	S1 <- stats::plogis( thetaM + extreme.item ) * 
				( 1 - stats::plogis( thetaM - extreme.item ) )
	likelihood <- likelihood * S1
	
# Revalpr( "mean(abs( like0 - likelihood) )")

	# likelihood adjustment
	like2 <- likelihood_adjustment_compute( likelihood, theta , thetaM , adjfac )
	

	#*** compute posterior given likelihood and empirical prior		
	if ( ! is.null( target.EAP.rel ) ){
		probs <- prob.theta	
		probsM <- matrix( prob.theta , nrow=N , ncol=TP , byrow=TRUE)

		tuning1 <- min_tuning
		tuning2 <- max_tuning
		
		
		EAP_rel1 <- likelihood_adjustment_tuning( likelihood , theta , thetaM , adjfac , 
						   tuningfac=tuning1 , probs=probs , normal_approx=normal_approx )
		EAP_rel2 <- likelihood_adjustment_tuning( likelihood , theta , thetaM , adjfac , 
						   tuningfac=tuning2 , probs=probs , normal_approx=normal_approx)	

		iter <- 0		
		change <- 1	
		while( ( iter < maxiter ) & ( change > conv ) ){									
			tuning0 <- ( tuning1 + tuning2 ) / 2									
			res1 <- likelihood_adjustment_tuning( likelihood , theta , thetaM , adjfac , 
							tuningfac=tuning0 , probs=probs , normal_approx=normal_approx)	
			EAP_rel0 <- res1$EAP.rel
			like2 <- res1$likelihood	
			if ( EAP_rel0 < target.EAP.rel ){
					tuning2 <- tuning0
						} else {
					tuning1 <- tuning0 
								}							
			iter <- iter + 1
			change <- abs( EAP_rel0 - target.EAP.rel )		 
			cat("Iteration ", iter , " | EAP reliability =" , round( EAP_rel0 , 4 ) , "\n")		
			flush.console()
		}	
			}
	res <- like2
	attr( res , "theta" ) <- matrix( theta , ncol=1)	
	attr(like0,"prob.theta") -> attr( res , "prob.theta")
	attr(like0,"G") -> attr(res , "G")	
	return(res)	
	}
####################################################################
	