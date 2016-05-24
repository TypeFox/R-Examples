##############################################
# MCMC estimation 2PNO model
mcmc.2pno <- function(dat , weights=NULL , burnin=500 , iter=1000 , N.sampvalues = 1000 ,
		progress.iter=50 , save.theta=FALSE ){
	s1 <- Sys.time()
	# data preparation
	dat0 <- dat
	dat <- as.matrix(dat)
	dat[ is.na(dat0) ] <- 0
	dat.resp <- 1-is.na(dat0)
	N <- nrow(dat)
	I <- ncol(dat)
	eps <- 10^(-10)
	# redefine weights
	if (! is.null(weights) ){
		weights <- N * weights / sum(weights )
					}	
	# set initial values
	a <- rep(1,I)
	b <- - stats::qnorm( (colMeans(dat0 , na.rm=TRUE) + .01 )/1.02 )
	# item parameters in matrix form
	aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
	bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)			
	theta <- qnorm( ( rowMeans( dat0,na.rm=TRUE ) + .01 ) / 1.02 )
	# define lower and upper thresholds 
	ZZ <- 1000
	threshlow <- -ZZ + ZZ*dat
	threshlow[ is.na(dat0) ] <- -ZZ
	threshupp <- ZZ*dat
	threshupp[ is.na(dat0) ] <- ZZ
	# saved values
	SV <- min( N.sampvalues  , iter - burnin )
    svindex <- round( seq( burnin , iter , len=SV )	 )
	a.chain <- matrix( NA , SV , I )
	b.chain <- matrix( NA , SV , I )
	theta.chain <- matrix( NA , SV , N )	
	deviance.chain <- rep(NA, SV)
	zz <- 0	
	#**********************
	# begin iterations
	for (ii in 1:iter){
		
		#****
		# draw latent data Z
		Z <- .draw.Z.2pl( aM , bM, theta , N , I , threshlow , threshupp )

		#***
		# draw latent traits theta
		res <- .draw.theta.2pl( aM , bM , N , I , Z )		
		theta <- res$theta

		#***
		# draw item parameters
		res <- .draw.itempars.2pl( theta , Z , I , N , weights)
		a <- res$a
		b <- res$b
		
		# item parameters in matrix form
		aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
		bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)		
				
		# save parameters
		if ( ii %in% svindex ){
			zz <- zz+1
			a.chain[ zz , ] <- a
			b.chain[ zz , ] <- b			
			theta.chain[ zz , ] <- theta			
			deviance.chain[zz] <- .mcmc.deviance.2pl( aM , bM , theta , dat , 
		     		dat.resp , 	weights , eps )
						}
		# print progress
		if ( ( ii %% progress.iter ) == 0 ){ 
			cat( "Iteration" , ii , " | " , paste(Sys.time()) , "\n")
			flush.console() }
						
				}
	##############################
	# output		
	# Information criteria
	ic <- .mcmc.ic.2pl( a.chain , b.chain , theta.chain , N , I , 
		dat , dat.resp , weights , eps ,
		deviance.chain )	
		
	# EAP reliability and person parameter estimates
	res <- .mcmc.person.2pno( theta.chain, weights )
	EAP.rel <- res$EAP.rel
	person <- res$person
		
	#-----
	# create MCMC object
	a <- a.chain 
	b <- b.chain
	theta <- theta.chain
	colnames(a) <- paste0("a[", 1:I , "]")
	colnames(b) <- paste0("b[", 1:I , "]")
	colnames(theta) <- paste0("theta[", 1:N , "]")
	mcmcobj <- cbind( deviance.chain , a , b )
	colnames(mcmcobj)[1] <- "deviance"
	if (save.theta){ mcmcobj <- cbind( mcmcobj , theta ) }
	class(mcmcobj) <- "mcmc"
	attr(mcmcobj, "mcpar") <- c( burnin+1 , burnin+SV , 1 )
	mcmcobj <- coda::as.mcmc.list( mcmcobj )
	
	#----
	# summary of the MCMC output
	summary.mcmcobj <- mcmc.list.descriptives( mcmcobj )
	# number of estimated parameters
#	np <- 2*I
	
	# time
	s2 <- Sys.time()
	time <- list( "start"=s1 , "end"=s2 , "timediff"= s2-s1 )
	#----
	# result list
	res <- list( "mcmcobj" = mcmcobj , "summary.mcmcobj" = summary.mcmcobj , 
			"ic"=ic , 
			"burnin"=burnin , "iter"=iter , 
			"a.chain"=a.chain , "b.chain" = b.chain , "theta.chain" = theta.chain ,
			"deviance.chain"=deviance.chain , 
			"EAP.rel" = EAP.rel , "person" = person , "dat" = dat0 , "weights" = weights , 
			"time" = time , "model" = "2pno" , "description"="2PNO Model" )
	class(res) <- "mcmc.sirt"
	return(res)				
	}