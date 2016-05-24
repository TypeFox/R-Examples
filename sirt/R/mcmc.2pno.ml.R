##############################################
# MCMC estimation 2PNO model
mcmc.2pno.ml <- function(dat , group , 
		link = "logit" , 
		est.b.M = "h" , est.b.Var = "n" , 
		est.a.M = "f" , est.a.Var = "n" , 
		burnin=500 , iter=1000 , N.sampvalues = 1000 ,
		progress.iter=50  , prior.sigma2 = c(1,.4 ) , 
		prior.sigma.b=c(1,1) , prior.sigma.a=c(1,1) ,
		prior.omega.b=c(1,1) , prior.omega.a=c(1,.4)  ,
		sigma.b.init=.3		){
	##############################
	# INPUT:	
	# describe input in help files
	##############################
	s1 <- Sys.time()
	# data preparation
	dat0 <- dat
	dat <- as.matrix(dat)
	dat[ is.na(dat0) ] <- 0
	dat.resp <- 1-is.na(dat0)
	N <- nrow(dat)
	I <- ncol(dat)
	eps <- 10^(-10)
	#***
	# groups
	groups <- sort( unique(group) )
	idgroup <- as.numeric( match( group , groups ) )
	G <- length(groups)
	groupsize <- as.numeric( rowsum( 1+0*group , idgroup )[,1] )
	
	#***
	# set initial values
	a <- rep(1,I)
	if (link=="logit"){ 
		b <- - stats::qnorm( (colMeans(dat0 , na.rm=TRUE) + .01 )/1.02 )
						}
	if (link=="normal"){b <- - colMeans( dat0 , na.rm=TRUE ) }
	eps.bG <- bG <- matrix( stats::rnorm( G*I , sd = sigma.b.init ) ,  G , I )	
	aG <- eps.aG <- 0*bG	
	sigma.res <- rep(1,I)
	if (link=="normal"){ 
		sigma.res <- apply( dat0 , 2 , stats::sd , na.rm=TRUE)
		Z <- dat
				}
	sigma.b <- rep( sigma.b.init , I )	
	sigma.a <- rep( .5 , I )
	# hyperparameters for items
	omega.b <- if ( est.b.M=="h"){ 2 } else { 1000 }
	mu.b <- 0
	omega.a <- if ( est.a.M=="h"){ 1 } else { 1000 }
	mu.a <- 1	
	# item parameters in matrix form
	aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
	bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE) + bG[ idgroup , ]
	aM.chainsum <- bM.chainsum <- 0*bM
	if (link=="logit"){	
		theta <- stats::qnorm( ( rowMeans( dat0,na.rm=TRUE ) + .01 ) / 1.02 )
						}
	if (link=="normal"){	
		theta <-  ( rowMeans( dat0,na.rm=TRUE ) + .01 ) / 1.02 
						}
	theta <- theta - mean( theta )					
	# theta level 2 part
	theta2 <- rep( stats::rnorm( G , sd = .2 ) , groupsize )	
	# inits Level 1 and Level 2 standard deviations
	sigma1 <- sqrt( .7 )
	sigma2 <- sqrt( .3 )
	# residual standard deviation (only applies in model with normal data)
	sigma.res <- apply( dat0 , 2  , stats::sd , na.rm=TRUE ) * .8
	#***
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
	sigma.res.chain <- sigma.a.chain <- sigma.b.chain <- matrix( NA , SV , I )
	theta.chain <- matrix( NA , SV , N )	
	theta2.chain <- matrix( NA , SV , G )
	# standard deviations
	sigma2.chain <- sigma1.chain <- rep( NA , SV )	# theta Level 1
	omega.a.chain <- mu.b.chain <- omega.b.chain <- rep(NA, SV )
	
	deviance.chain <- rep(NA, SV)
	zz <- 0	
	#**********************
	# begin iterations
	for (ii in 1:iter){
		
		#****
		# draw latent data Z
		if (link=="logit"){
			Z <- .draw.Z.2pno.ml( aM , bM, theta , N , I , threshlow , threshupp )
						}

		#***
		# draw latent traits theta (total)
		theta <- .draw.theta.2pno.ml( aM , b , bM , N , I , Z ,
				sigma1 , sigma2 , sigma.res , link , theta2 , idgroup  )		

		#***
		# draw latent group means (level 2 effects)
		theta2 <- .draw.theta2.2pno.ml( theta , idgroup , groupsize , 
						sigma1 , sigma2 , G )	
		
		#***
		# draw level 1 and level 2 theta variances
		res <- .draw.sigma12.2pno.ml( theta , theta2 , idgroup , N ,
				G , prior.sigma2 )		
		sigma1 <- res$sigma1
		sigma2 <- res$sigma2

		
		#***
		# draw item groupwise parameters b (est.b.Var = "j" or est.b.Var = "i")
		if ( est.b.Var %in% c("i","j")){

			# sampling of b parameters		
			b <- .mcmc.est.b.2pno.ml.v2( N ,Z , aM , theta , idgroup , groupsize , 
										b , bG , G , I , sigma.b , omega.b , 
										mu.b , sigma.res , link )			
			# sampling of groupwise b item parameters
			bG <- .mcmc.est.b.group.2pno.ml( Z , aM , theta , idgroup , groupsize , 
					b , G , I , sigma.b , sigma.res , link , N  )
			# sampling of sigma.b
			sigma.b <- .mcmc.sigma.b.2pno.ml( bG , mu.b , omega.b , G , I , 
				est.b.Var , prior.sigma.b , sigma.b )			
			# define item parameters in a matrix								
			bM <- matrix(b,N,I,byrow=TRUE) + bG[ idgroup , ]

				}				
		#***
		# draw item difficulties (est.b.Var = "n")
		# single level case
		if ( est.b.Var == "n"){
			b <- .draw.est.b.sl(Z , aM , theta , N , I , omega.b , mu.b ,
				sigma.res , link )
		    bM <- matrix( b , nrow=N , ncol=I , byrow=TRUE)					
						}
						
		# draw hyperparameters of b (independent of est.b.Var)
		res <- .draw.est.b.hyperpars( b, mu.b , omega.b , I , 
			      prior.omega.b , est.b.M  )
		mu.b <- res$mu.b
		omega.b <- res$omega.b
		
		##################################
		# sampling of a parameters
		if ( est.a.M != "f" ){
			#****
			# draw a item parameters group wise
			if ( est.a.Var %in% c("i","j")){

				# sampling of a parameters
				a <- .mcmc.a.est.a.2pno.ml( Z , bM , aG , idgroup , theta , 
						mu.a , omega.a , I , link , sigma.res  )				
				
				# sampling of groupwise a parameters
				aG <- .mcmc.est.aG.2pno.ml.v2( Z , bM , theta , idgroup , G , I ,
						a , sigma.a , N , link , sigma.res )
						
				# sampling of standard deviations of random item 
				# discrimination effects
				sigma.a <- .mcmc.a.grouphier.2pno.ml( aG , mu.a , G , omega.a , I ,
							prior.sigma.a , est.a.Var , sigma.a )
				# calculate total item discrimination
				aM <- matrix( a , N , I , byrow=TRUE ) + aG[ idgroup , ]			
											}		
			#***
			# draw a item parameters (single level estimation)
			if ( est.a.Var =="n"){
				a <- .draw.est.a.sl( Z , bM , theta , mu.a , omega.a , I ,
						sigma.res , link)		
				aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE)
						}
						
			# draw hyperparameters of a
			res <- .draw.est.a.hyperpars( a, mu.a , omega.a , I , 
					prior.omega.a , est.a.M  )
			omega.a <- res$omega.a
							}

		# draw residual standard deviations for link="normal"	
		if (link=="normal"){
			sigma.res <- .draw.sigma.res.2pno.ml( Z , aM , bM , theta , N , I )
					}							
							
		# save parameters
		if ( ii %in% svindex ){
			zz <- zz+1
			a.chain[ zz , ] <- a
			b.chain[ zz , ] <- b			
			theta.chain[ zz , ] <- theta			
			theta2.chain[ zz , ] <- theta2	
			if (link=="logit"){			
				deviance.zz <- .mcmc.deviance.2pl( aM , bM , theta , dat , 
						dat.resp , 	weights=NULL , eps )
							}
			if (link=="normal"){
				deviance.zz <- .mcmc.deviance.normallink.2pno.ml( aM , bM , 
						theta , N , I , dat , dat.resp , sigma.res )
								}	
			deviance.chain[zz] <- deviance.zz
			sigma1.chain[zz] <- sigma1
			sigma2.chain[zz] <- sigma2
			sigma.b.chain[zz,] <- sigma.b
			sigma.a.chain[zz,] <- sigma.a
			mu.b.chain[zz]	<- mu.b
			omega.b.chain[zz] <- omega.b
			omega.a.chain[zz] <- omega.a	
			if (link=="normal"){
				sigma.res.chain[zz,] <- sigma.res
					}
			# sum b and a parameters
			bM.chainsum <- bM.chainsum + bM
			aM.chainsum <- aM.chainsum + aM
						}
		# print progress
		if ( ( ii %% progress.iter ) == 0 ){ 
			cat( "Iteration" , ii , " | " , paste(Sys.time()) , "\n")
			flush.console() }
						
				}
	##############################
	# output		
	# Information criteria
	ic <- .mcmc.ic.2pno.ml( aM.chainsum , bM.chainsum , theta.chain , N , I , 
				dat , dat.resp ,eps , deviance.chain , groupsize ,
				sigma.res.chain , link  )
		
	# EAP reliability and person parameter estimates (okay!)
	res <- .mcmc.person.2pno.ml( theta.chain, weights=NULL )
	EAP.rel <- res$EAP.rel
	person <- res$person
		
	#-----
	# create MCMC object
	mcmcobj <- .mcmc.mcmclist.2pno.ml( a.chain , b.chain , deviance.chain ,
		sigma1.chain , sigma2.chain , I , SV , burnin , iter ,
		est.b.M , mu.b.chain , omega.b.chain , sigma.b.chain , est.b.Var ,
		est.a.M , est.a.Var , omega.a.chain , sigma.a.chain ,
		sigma.res.chain , link		)
	
	#----
	# summary of the MCMC output
	summary.mcmcobj <- mcmc.list.descriptives( mcmcobj )
	
	# time
	s2 <- Sys.time()
	time <- list( "start"=s1 , "end"=s2 , "timediff"= s2-s1 )
	# description
	if ( link=="logit"){ description <- "2PNO Multilevel Model" }
	if ( link=="normal"){ description <- "Normal Multilevel Model" }	
	#----
	# result list
	res <- list( "mcmcobj" = mcmcobj , "summary.mcmcobj" = summary.mcmcobj , 
			"ic"=ic , 
			"burnin"=burnin , "iter"=iter , 
			"theta.chain" = theta.chain , "theta2.chain" = theta2.chain , 
			"deviance.chain"=deviance.chain , 
			"EAP.rel" = EAP.rel , "person" = person , "dat" = dat0 , 
			"time" = time , "model" = "2pno.ml" , "description"= description )
	class(res) <- "mcmc.sirt"
	return(res)				
	}