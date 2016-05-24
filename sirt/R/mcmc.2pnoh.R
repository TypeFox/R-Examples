##############################################
# MCMC estimation 2PNO model
mcmc.2pnoh <- function(dat , itemgroups , prob.mastery = c(.5 , .8 ) ,
		weights=NULL , burnin=500 , iter=1000 , N.sampvalues = 1000 ,
		progress.iter=50 , prior.variance=c(1,1) , save.theta=FALSE ){
	s1 <- Sys.time()
	# data preparation
	dat0 <- dat
	dat <- as.matrix(dat)
	dat[ is.na(dat0) ] <- 0
	dat.resp <- 1-is.na(dat0)
	N <- nrow(dat)
	I <- ncol(dat)
	eps <- 10^(-10)
	# itemgroups
	itemgroups.unique <- sort( unique( itemgroups ) )
	K <- length(itemgroups.unique)
	itemgroup <- match( itemgroups , itemgroups.unique )
	Ik <- stats::aggregate( 1 + 0*1:I , list(itemgroup) , sum )[,2]
	# redefine weights
	if (! is.null(weights) ){
		weights <- N * weights / sum(weights )
					}	
	# set initial values
	a <- rep(1,I)
	b <- - stats::qnorm( (colMeans(dat0 , na.rm=TRUE) + .01 )/1.02 )
	xi <- rep(0,K)
	omega <- rep(1,K)
	sig <- 1.5	# SD of item difficulties
	nu <- .30	# SD of item discriminations
	# item parameters in matrix form
	aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
	bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)			
	theta <- stats::qnorm( ( rowMeans( dat0,na.rm=TRUE ) + .01 ) / 1.02 )
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
	nu.chain <- sig.chain <- deviance.chain <- rep(NA, SV)
	transition.chain <- nonmastery.chain <- mastery.chain <- 
			xi.chain <- omega.chain <- matrix(NA , SV , K )
	M.beta.chain <- M.alpha.chain <- matrix(NA,SV,K)
	zz <- 0	
	#**********************
	# begin iterations
	for (ii in 1:iter){
		
		#****
		# draw latent data Z
		Z <- .draw.Z.2pnoh( aM , bM, theta , N , I , threshlow , threshupp )
		
		#***
		# draw latent traits theta
		res <- .draw.theta.2pnoh( aM , bM , N , I , Z )		
		theta <- res$theta
			
		#***
		# draw item parameters alpha, beta, xi and omega
		res <- .draw.itempars.2pnoh( theta , Z , I , N , a , b , 
				xi , omega , sig , nu , itemgroup , K , Ik , weights)
		a <- res$a
		b <- res$b
		xi <- res$xi
		omega <- res$omega
		
		# draw item variance
		res <- .draw.itemvariances.2pnoh(a,b,I , itemgroup , xi , omega , prior.variance)
		sig <- res$sig
		nu <- res$nu

		# item parameters in matrix form
		aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
		bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)		
					
		# save parameters
		if ( ii %in% svindex ){
			zz <- zz+1
			a.chain[ zz , ] <- a
			b.chain[ zz , ] <- b
			xi.chain[zz,] <- xi
			omega.chain[zz,] <- omega
			theta.chain[ zz , ] <- theta		
			# mean alpha and beta
#			M.alpha.chain[zz,] <- aggregate( a , list(itemgroup) , mean)[,2]
			M.alpha.chain[zz,] <- rowsum( a , itemgroup )[,1] / Ik
			M.beta.chain[zz,] <- rowsum( b , itemgroup )[,1] / Ik
			# calculate deviance
			deviance.chain[zz] <- .mcmc.deviance.2pl( aM , bM , theta , dat , 
		     		dat.resp , 	weights , eps )
			# calculate mastery probabilities
			res <- .mcmc.mastery.2pnoh( xi , omega , N , K , weights , theta , prob.mastery)
			nonmastery <- res$nonmastery 
			transition <- res$transition
			mastery <- res$mastery					
			nonmastery.chain[zz,] <- res$nonmastery
			transition.chain[zz,] <- res$transition
			mastery.chain[zz,] <- res$mastery			
			nu.chain[zz] <- nu
			sig.chain[zz] <- sig
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
	# MCMC object
	a <- a.chain 
	b <- b.chain
	d <- a*b
	theta <- theta.chain
	colnames(a) <- paste0("alpha[", 1:I , "]")
	colnames(b) <- paste0("beta[", 1:I , "]")
	colnames(d) <- paste0("d[", 1:I , "]")
	colnames(theta) <- paste0("theta[", 1:N , "]")
	mcmcobj <- cbind( deviance.chain , a , b , d )
	colnames(mcmcobj)[1] <- "deviance"	
	colnames(xi.chain) <- paste0("xi[", 1:K , "]")	
	colnames(omega.chain) <- paste0("omega[", 1:K , "]")
	colnames(M.alpha.chain) <- paste0("M.alpha[", 1:K , "]")	
	colnames(M.beta.chain) <- paste0("M.beta[", 1:K , "]")	
	tau.chain <- omega.chain * xi.chain	
	colnames(tau.chain) <- paste0("tau[", 1:K , "]")	
	mcmcobj <- cbind( mcmcobj , xi.chain , omega.chain , 
		M.beta.chain , M.alpha.chain , tau.chain )
	dfr <- cbind( "sig" = sig.chain , "nu"=nu.chain )
	mcmcobj <- cbind( mcmcobj , dfr )
	colnames(nonmastery.chain) <- paste0("nonmastery[", 1:K , "]")	
	colnames(transition.chain) <- paste0("transition[", 1:K , "]")	
	colnames(mastery.chain) <- paste0("mastery[", 1:K , "]")	
	mcmcobj <- cbind( mcmcobj , nonmastery.chain , transition.chain , mastery.chain )
	if (save.theta){ mcmcobj <- cbind( mcmcobj , theta ) }
	class(mcmcobj) <- "mcmc"
	attr(mcmcobj, "mcpar") <- c( burnin+1 , burnin+SV , 1 )
	mcmcobj <- coda::as.mcmc.list( mcmcobj )
	
	#----
	# summary of the MCMC output
	summary.mcmcobj <- mcmc.list.descriptives( mcmcobj )
	# number of estimated parameters
#	np <- 2*I + 2*K + 2
	
	s2 <- Sys.time()
	time <- list( "start"=s1 , "end"=s2 , "timediff"= s2-s1 )
	#----
	# result list
	res <- list( "mcmcobj" = mcmcobj , "summary.mcmcobj" = summary.mcmcobj , 
			"ic"=ic , 
			"burnin"=burnin , "iter"=iter , 
			"alpha.chain"=a.chain , "beta.chain" = b.chain , 
			"xi.chain" = xi.chain , "omega.chain" = omega.chain ,
			"sig.chain" = sig.chain , "nu.chain" = nu.chain , 
			"theta.chain" = theta.chain ,
			"deviance.chain"=deviance.chain , 
			"EAP.rel" = EAP.rel , "person" = person , 
			"dat" = dat0 , "weights" = weights , 
			"time" = time , "model" = "2pnoh" ,
			"description"="2PNO Hierarchical IRT Model for Criterion-Referenced Measurement")
	class(res) <- "mcmc.sirt"
	return(res)				
	}