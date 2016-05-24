
#########################################################
# Functional Unidimensional Model (Ip et al., 2013)
f1d.irt <- function( dat=NULL , nnormal=1000 , nfactors=3 ,
			A=NULL , intercept=NULL , mu=NULL , Sigma = NULL , maxiter=100 ,
			conv=10^(-5) , progress=TRUE ){
	if ( ! is.null(dat) ){		
		# estimate tetrachoric correlation matrix
		if (progress){
			cat("*** Estimate tetrachoric correlation\n")
					}
		tetra <- res <- tetrachoric2(dat , progress=progress)
		# estimate factor analysis
		fac1 <- psych::fa( r=res$rho , nfactors=nfactors , rotate="none" )
		fac0 <- psych::fa( r=res$rho , nfactors=1 , rotate="none" )
		# extract standardized loadings
		A_stand <- as.matrix( fac1$loadings )
		# calculate communality
		h2 <- rowSums( A_stand^2 )
		# unstandardized loadings
		A <- A_stand / sqrt( 1-h2 )
		# intercepts
		intercept <- - res$tau / sqrt( 1 - h2 )
		names.dat <- colnames(dat)
					} else {
		a0 <- NA
		d0 <- NA
		A[ is.na(A) ] <- 0
		names.dat <- names(intercept)
		tetra <- NULL
				}
						
	#***************************************
	# approximation of normal distribution using quasi Monte Carlo integration nodes
    theta <- qmc.nodes( nnormal , nfactors )
	if ( is.null(mu) ){
		mu <- rep(0,nfactors)
						}	
	if ( is.null(Sigma) ){
		Sigma <- diag(1,nfactors)
						}
	
	wgt_theta <- mvtnorm::dmvnorm(x=theta, mean= mu, sigma= Sigma )
	wgt_theta <- wgt_theta / sum( wgt_theta )

	I <- length(intercept)
	TP <- nrow(theta)
	W1 <- diag(wgt_theta )

	D <- ncol(A)
	# a_i ' theta_p
	Zpi <- matrix( 0 , TP , I )
	for (dd in 1:D){
		# dd <- 1
		Zpi <- Zpi + theta[,dd] * matrix( A[,dd] , TP , I , byrow=TRUE )
				  }
	# Z_pi = a_i theta_p + d_i
	Zpi <- Zpi + matrix( intercept , TP , I , byrow=TRUE )

	# starting values for a_i and d_i
	di <- intercept
	diast <- di  # diast = d_i*
	aiast <- rep(.7,I)

	thetaast <- rep(0,TP)
	iter <- 0
	parchange <- 1

	#*****************************************
	# begin algorithm
	while( ( iter < maxiter ) & ( parchange > conv ) ){
	
		thetaast0 <- thetaast
		aiast0 <- aiast
		diast0 <- diast

		# (1) update theta_p using ( a_i , theta_p , d_i , a_i* , d_i* )
		Ypi <- Zpi - matrix( diast , TP , I , byrow=TRUE)
		aiastM <- matrix( aiast , TP , I , byrow=TRUE )
		thetaast <- rowSums( Ypi * aiastM ) / rowSums( aiastM^2 )
		wM <- stats::weighted.mean( thetaast , wgt_theta )
		sdM <- sqrt( sum( ( thetaast - wM )^2 * wgt_theta ) )
		thetaast <- ( thetaast - wM ) / sdM

		X <- cbind( 1 , thetaast )
		c1 <- solve( crossprod( X , W1 ) %*% X , crossprod( X , W1 ) %*% Zpi )
		diast <- c1[1,]
		aiast <- c1[2,]            

		# compute approximation error
		errpi <- Zpi - thetaast * matrix(aiast,TP , I , byrow=TRUE )  - 
				 matrix(diast,TP , I , byrow=TRUE )
		approx.error <- sum( errpi^2 * wgt_theta ) / I
		# parameter change
		parchange <- max( abs( c(diast - diast0,aiast-aiast0,thetaast-thetaast0) ))

		iter <- iter + 1
        if (progress){
			cat( paste0( "Iteration " , iter , 
				" | Approximation error = " , round( approx.error , 5 ) ,
				" | Max. parameter change = " , round( parchange , 5) ,
				 "\n") )
			utils::flush.console()
					}
		}
	#**************************************************
	

	if ( ! is.null(dat) ){
		# unstandardized loadings 1 factor model
		A0_stand <- fac0$loadings
		a0 <- A0_stand[,1] / sqrt( 1 - A0_stand[,1]^2 )
		d0 <- - res$tau / sqrt( 1 - A0_stand[,1]^2 )
					}
	if ( is.null(dat) ){
		a0 <- NULL
		d0 <- NULL
				}

	item <- data.frame(  "item" = names.dat )
	item$ai.ast <- aiast 
	item$ai0 <- a0 
	item$di.ast <- diast  
	item$di0 <- d0 

	person <- data.frame( "theta.ast" = thetaast ,
			"wgt" = wgt_theta )

	res <- list( "item" = item , "person"=person , 
				"A"=A , "intercept"=intercept ,
				"dat"=dat , "tetra" = tetra )
	return(res)
		}
#**************************************************