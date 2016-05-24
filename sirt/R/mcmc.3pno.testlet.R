######################################################
# MCMC estimation 2PNO model
mcmc.3pno.testlet <- function(dat , testlets=rep(NA , ncol(dat)) , 
		weights=NULL , est.slope =TRUE , est.guess = TRUE , guess.prior = NULL , 
		testlet.variance.prior = c( 1 , .2 ) , 
		burnin=500 , iter=1000 , N.sampvalues = 1000 ,
		progress.iter=50 , save.theta=FALSE ){
	s1 <- Sys.time()
	param <- 1
#	if (param>2){	
#		cat("Implementation of parametrization 3 is broken and invalid\n")
#		stop("Use only parameterizations 1 or 2!\n")
#			}	
	# data preparation
	dat0 <- dat
	dat <- as.matrix(dat)
	dat[ is.na(dat0) ] <- 0
	dat.resp <- 1-is.na(dat0)
	N <- nrow(dat)
	I <- ncol(dat)
	eps <- 10^(-10)
	# testlet groups
	testlets.u <- setdiff( unique( testlets ) , c(NA) )
	TT <- length(testlets.u)
	testletgroups <- match( testlets , testlets.u )
    testletgroups[ is.na(testletgroups ) ] <- TT+1
	N.items <- table( testletgroups )	
		
	# redefine weights
	if (! is.null(weights) ){
		weights <- N * weights / sum(weights )
					}	
	# set initial values
#	a.testlet <- a <- rep(.6,I)	
	a.testlet <- a <- rep(1,I)	
	b <- stats::qnorm( (colMeans(dat0 , na.rm=TRUE) + .01 )/1.02 )
	# guessing parameters
	if ( is.null(guess.prior) ){
		guess.prior <- matrix( c(1,1) , nrow=I , ncol=2 , byrow=TRUE )
		guess <- rep( .4 , I )
								} else {
		if ( is.vector( guess.prior) ){
			guess.prior <- matrix( guess.prior , nrow=I , ncol=2 , byrow=TRUE )
							}		
		guess <- guess.prior[,1] / rowSums( guess.prior )
					}				
	if ( ! est.guess ){ guess <- rep(0,I) }					
	# item parameters in matrix form
	a.testletM <- aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
	bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)			
	# theta effects
	theta <- stats::qnorm( ( rowMeans( dat0,na.rm=TRUE ) + .01 ) / 1.02 )
	theta <- theta - mean(theta)
	mu <- 0
	sigma <- 1
	# testlet effects
	gamma.testlet <- matrix( rnorm( N*(TT+1) , sd = .4 ) , N , TT+1 )
	gamma.testlet[,TT+1] <- 0	
	# testlet variances
	sigma.testlet <- rep( .4 , TT )
	if (param==3){ sigma.testlet <- rep(1,TT) }
	# define lower and upper thresholds 
	ZZ <- 1000
	threshlow <- -ZZ + ZZ*dat
	threshlow[ is.na(dat0) ] <- -ZZ
	threshupp <- ZZ*dat
	threshupp[ is.na(dat0) ] <- ZZ
	W <- dat
	# saved values
	SV <- min( N.sampvalues  , iter - burnin )
    svindex <- round( seq( burnin , iter , len=SV )	 )
	a.chain <- matrix( NA , SV , I )
	a.testlet.chain <- matrix(NA,SV,I)
	sigma.testlet.chain <- matrix(NA , SV , TT )
	c.chain <- b.chain <- matrix( NA , SV , I )
	theta.chain <- matrix( NA , SV , N )	
	sigma.chain <- deviance.chain <- rep(NA, SV)
	gamma.testlet.chain <- matrix(NA, SV , N*(TT+1) )
	zz <- 0	
		
	
	#************************************************	
	#************************************************
	# begin iterations
	for (ii in 1:iter){
		
		#************
		# param = 1
		# tau_{ni} = a_i \theta_n + \gamma_{nt} + b_i
		# param = 2
		# tau_{ni} = a_i \theta_n + a_i \gamma_{nt} + b_i		

		# define some intermediate parameters
		gamma.testletM <- gamma.testlet[ , testletgroups ]	
		if (param==1){ tau.ni <- aM * theta + gamma.testletM + bM }
		if (param==2){ tau.ni <- aM * theta + aM * gamma.testletM + bM }		
		if (param==3){ tau.ni <- aM * theta + a.testletM * gamma.testletM + bM }		
		
		#***
		# draw latent data W
		if (est.guess){
			W <- .draw.W.3pno.testlet( aM , bM, theta , gamma.testlet ,
					N , I , threshlow , threshupp , testletgroups , param , dat ,
					guess , tau.ni)
					}
		#****
		# draw latent data Z
		Z <- .draw.Z.3pno.testlet( aM , bM, theta , gamma.testlet ,
				N , I , threshlow , threshupp , testletgroups , param , W ,
				dat.resp , tau.ni)
			
		#***
		# draw latent traits theta
		theta <- .draw.theta.3pno.testlet( aM , bM , N , I , Z , param , gamma.testletM ,
			     sigma ,  a.testletM)		 
			 
		# draw theta SD
		if ( ! est.slope ){ 
			sigma <- .draw.theta.variance.3pno.testlet( theta , weights , N )
						   }				 
		#***
		# draw testlet effects
		if (TT>0){
			gamma.testlet <- .draw.gamma.3pno.testlet( aM , bM , N , I , Z , param ,
				theta , testletgroups , TT , N.items , gamma.testlet , 
				sigma.testlet , a.testletM)		
				}			
		#***
		# draw item parameters a and b
		if ( est.slope ){
			res <- .draw.itempars.3pno.testlet( theta , Z , I , N , weights ,
					gamma.testlet , testletgroups , param  , TT , a.testletM )
			a <- res$a ; b <- res$b		
						}						
		# draw only b item parameters
		if ( ! est.slope ){
			b <- .draw.est.b.3pno.testlet( aM , bM , N , I , Z , param ,
					gamma.testletM , sigma , weights , theta )
						}
									
		#***
		# draw guessing parameter
		if (est.guess){
			guess <- .draw.guess.3pno.testlet( guess.prior , W , 
					dat , dat.resp , weights , I )
					   }
		# draw testlet standard deviations
		if ( (TT>0) & (param!=3) ){
			sigma.testlet <- .draw.testlet.variance.3pno.testlet( gamma.testlet , N ,
					   sigma.testlet , testlet.variance.prior , weights , TT)		
				  }
		# draw testlet slope
		if (param==3){	 		
			a.testlet <- .draw.est.a.testlet.3pno.testlet( aM , bM , N , I , Z , param ,
					gamma.testletM , sigma , weights , theta , testletgroups ,
					gamma.testlet , TT	)		
						}
				  
		# define item parameters in matrix form
		aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
		bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)		
		if (param==3){ a.testletM <- matrix( a.testlet , nrow=N , ncol=I , byrow=TRUE) }
	
		
		# save parameters
		if ( ii %in% svindex ){
			zz <- zz+1
			a.chain[ zz , ] <- a					
			b.chain[ zz , ] <- b			
			c.chain[ zz , ] <- guess
			if (param==3){
				a.testlet.chain[zz,] <- a.testlet
						}
			sigma.testlet.chain[zz,] <- sigma.testlet
			theta.chain[ zz , ] <- theta
			sigma.chain[ zz ] <- sigma
			gamma.testlet.chain[zz , ] <- matrix(gamma.testlet , ncol=N*(TT+1) , nrow=1 )
			deviance.chain[zz] <- .mcmc.deviance.3pno.testlet( aM , bM , 
					a.testletM , theta , guess ,	gamma.testlet , testletgroups , 
					dat , dat.resp , weights , eps , param )
						}
		# print progress
		if ( ( ii %% progress.iter ) == 0 ){ 
			cat( "Iteration" , ii , " | " , paste(Sys.time()) , "\n")
			flush.console() }					
				}   # end MCMC iterations
	##############################
	# output	

	# EAP reliability and person parameter estimates
	res <- .mcmc.person.3pno.testlet( theta.chain , weights ,
				gamma.testlet.chain )
	EAP.rel <- res$EAP.rel
	person <- res$person
	
	#----
	# calculate information criteria	
	ic <- .mcmc.ic.3pno.testlet( a.chain , b.chain , a.testlet.chain , N , I , 
		theta.chain , c.chain , gamma.testlet.chain , TT , 
		testletgroups , dat , dat.resp , weights , eps , param ,
		deviance.chain )	
	
	#-----
	# create MCMC object
	res <- .mcmc.list.3pno.testlet( a.chain , b.chain , a.testlet.chain ,
				I , deviance.chain , est.slope , c.chain , sigma.chain ,
				est.guess , theta.chain , sigma.testlet.chain , TT ,
				burnin , SV , save.theta , testletgroups , param )
	mcmcobj <- res$mcmcobj
	theta <- res$theta
	
	#----
	# summary of the MCMC output
	summary.mcmcobj <- mcmc.list.descriptives( mcmcobj )
	
	#****
	# description
	description <- .mcmc.description.3pno.testlet( TT , param , est.guess , est.slope )

	#***
	# time
	s2 <- Sys.time()
	time <- list( "start"=s1 , "end"=s2 , "timediff"= s2-s1 )
	#----
	# result list
	res <- list( "mcmcobj" = mcmcobj , "summary.mcmcobj" = summary.mcmcobj , 
			"ic"=ic , 
			"burnin"=burnin , "iter"=iter , 
			"theta.chain" = theta.chain ,
			"deviance.chain"=deviance.chain , 
			"EAP.rel" = EAP.rel , "person" = person , "dat" = dat0 , "weights" = weights , 
			"time" = time , "model" = "3pno.testlet" , "description"=description , 
			"TT"=TT  )
	class(res) <- "mcmc.sirt"
	return(res)				
	}