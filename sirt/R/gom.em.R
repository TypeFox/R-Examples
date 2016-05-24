
#################################################################
# gom EM algorithm
gom.em <- function( dat , K=NULL , problevels=NULL , model="GOM" , 
	theta0.k= seq(-5,5,len=15) , xsi0.k= exp(seq(-6,3 ,len=15)) , 
	max.increment=.3 , numdiff.parm=.001 , maxdevchange= 10^(-5) ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001 ,
	progress=TRUE ){
	#..........................................................
	s1 <- Sys.time()
	e1 <- environment() 
	s1 <- Sys.time()
	if ( model=="GOMRasch"){
		K <- length(theta0.k)
		max.increment <- 1
					}
#	res <- .gom.proc( dat , K =K ) 
#	.sirt.attach.environment( res , envir=e1 )
    dat0 <- dat 
    dat.resp <- 1-is.na(dat)
    dat[ is.na(dat) ] <- 0
    N <- nrow(dat)
    I <- ncol(dat)
	dat2 <- as.matrix(dat)
	dat2.resp <- as.matrix(dat.resp)
	# indicator matrix
	dat2.ind0 <- dat2.resp * 1*(dat2==0)
	dat2.ind1 <- dat2.resp * 1*(dat2==1)
	dat2.ind <- as.matrix( cbind( dat2.ind0 , dat2.ind1 ) )
	mu <- Sigma <- b <- NULL
	# design matrices
	if (model=="GOMRasch"){ 
		b <- - stats::qlogis( colMeans(dat, na.rm=T)  )
		theta.kM <- as.matrix( expand.grid( theta0.k , xsi0.k ))
		TP <- nrow(theta.kM)
		m1 <- exp( - ( theta.kM[,1] - matrix( theta0.k , TP , K , byrow=T ) )^2 / ( 2*theta.kM[,2]^2 ) )
		theta.k <- m1 / rowSums( m1 )		
		#***
		mu <- c(0, .7)
		Sigma <- as.matrix(diag(c(1 , 1 )))
		pi.k <- mvtnorm::dmvnorm( theta.kM , mu , Sigma )
		pi.k <- pi.k / sum( pi.k )
		lambda <- stats::plogis( - outer( b , theta0.k , "-" ) )
		# design matrix skill space
#		Z <- cbind( 1 , theta.kM[,1] , theta.kM[,1]^2 , theta.kM[,1]^3 ,
#				 theta.kM[,2] , theta.kM[,2]^2 , theta.kM[,2]^3 ,
#				 theta.kM[,1]*theta.kM[,2] , theta.kM[,1]^2*theta.kM[,2])	
				}
	
	if (model=="GOM"){
		theta.k <- .gom.calc.theta(K, problevels)
		TP <- nrow(theta.k)		
		pi.k <- as.vector(rep(1/TP , TP ))
		lambda <- matrix( .75*seq( 1/(2*K) , 1 , 1/K) , I , K , byrow=T ) 
		theta0.k <- NULL
					}	
	# inits
	se.b <- se.lambda <- NULL
	n.ik <- array( 0 , dim=c(TP,I , 2) )
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	disp <- "...........................................................\n"	
	
	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter )
						){
		if (progress){				
			cat(disp)	
			cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
					}
		# previous values
		dev0 <- dev
		pi.k0 <- pi.k
		lambda0 <- lambda
		b0 <- b
		# calculate probabilities
		res <- .gom.calcprobs( lambda , theta.k , b , theta0.k )   
		probs <- res$probs
		probs <- problong2probarray( probs , I , TP )
		
		# calculate counts
		probsM <- matrix( base::aperm( probs , c(2,1,3) ) , nrow=I*2 , ncol=TP )	
		res1 <- calcpost( dat2 , dat2.resp , probsM , dat2.ind , pi.k , K=1 )
		f.yi.qk <- res1$fyiqk
		f.qk.yi <- res1$f.qk.yi
		pi.k <- res1$pi.k
		n.ik <- array( res1$n.ik , dim=dim(n.ik) )		
		N.ik <- res1$N.ik

		# maximize lambda
		if (model=="GOM"){
			res <- .gom.est.lambda( lambda , I , K , n.ik , 
					   numdiff.parm , max.increment=max.increment, theta.k , msteps ,
					   mstepconv , eps = .001 , progress=progress )
			lambda <- res$lambda
			se.lambda <- res$se.lambda
			max.increment <- max( abs(lambda-lambda0))/1.2	
							}
		if (model=="GOMRasch"){
			res <- .gom.est.b( lambda , I , K , n.ik , b , theta0.k , 
					numdiff.parm=.001 , max.increment=max.increment,theta.k , msteps ,
					mstepconv , eps = .001 , progress=progress)
			b <- res$b
			se.b <- res$se.b
			lambda <- t( stats::plogis( outer( theta0.k , b , "-" )	) )
			max.increment <- max( abs(b-b0))/1.2	
							}
		flush.console()		
		
		#**
		# calculate deviance
#		pi.k <- matrix( pi.k , ncol=1 )	
		ll <- sum( log( rowSums( f.yi.qk * outer( rep(1,nrow(f.yi.qk)) , pi.k ) ) ) )
#		ll <- sum( log( rowSums( f.yi.qk * outer( rep(1,nrow(f.yi.qk)) , pi.k[,1] ) ) ) )
		dev <- -2*ll
		# convergence criteria
		conv <- max( abs(lambda-lambda0)) 
		iter <- iter+1
		devchange <- abs( ( dev - dev0 )/dev0 )	
		#****
		# print progress
        if (progress){		
			cat( paste( "   Deviance = "  , round( dev , 4 ) , 
				if (iter > 1 ){ " | Deviance change = " } else {""} ,
				if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,"\n",sep="") )
			cat( paste( "    Maximum lambda parameter change = " , 
					paste( round(max(abs(lambda-lambda0)) ,6) , collapse=" " ) , "\n" , sep=""))
			cat( paste( "    Maximum distribution parameter change = " , 
					paste( round(max(abs(pi.k-pi.k0)) ,6) , collapse=" " ) , "\n" , sep=""))
			if (model=="GOMRasch"){
				cat( paste( "    Maximum b parameter change = " , 
					paste( round(max(abs(b-b0)) ,6) , collapse=" " ) , "\n" , sep=""))
						}				
					}				
				}
	# *********
	# arrange OUTPUT
	#---
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$np.item <- I*K
	if (model=="GOMRasch"){ ic$np.item <- I }	
	# trait matrix
	ic$np.trait <- TP - 1
	if (model=="GOMRasch"){ ic$np.trait <- 3 + 1 }
	ic$np <- ic$np.item + ic$np.trait
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	
	#---
	# person parameters
	#---
	# item
	item <- data.frame("item" = colnames(dat))
	item$N <- colSums( dat2.resp )
	item$p <- colMeans( dat2 , na.rm=TRUE)
	item$b <- b
	if (model != "GOMRaschxxx"){
		for (kk in 1:K){
			item[,paste0("lam.Cl",kk)] <- lambda[,kk] 
					}
						}
	obji <- item
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
    if (progress){		
		cat("*********************************\n")
		cat("Item Parameters\n")
		print( obji )		
			}
	EAP.rel <- NULL
	person <- NULL
	if ( model=="GOMRasch"){			
		res <- .gom.est.covariance( f.qk.yi , Sigma , theta.kM , N  )
		.sirt.attach.environment( res , envir=e1 )				
		#--- distribution parameters
		c1 <- cov2cor(Sigma)
		if (progress){
			cat("*********************************\n")
			cat("Trait Distribution (Location, Variability)\n")		
			cat( " Means: " , round( mu , 3 ) , "\n")				
			cat( " Standard deviations: " , round( sqrt(diag(Sigma)) , 3 ) , "\n")
			cat( " Correlation " , round( c1[lower.tri(c1)] , 3 ) , "\n")		
			flush.console()			
					}
		# person parameters
		pers <- .smirt.person.parameters( data=dat2 , D=2 , theta.k=theta.kM ,
			p.xi.aj=f.yi.qk , p.aj.xi=f.qk.yi , weights=rep(1,N) )	
		person <- pers$person
		EAP.rel <- pers$EAP.rel 
		if (progress){
			cat("*********************************\n")
			cat("EAP Reliability = " , round(EAP.rel,3) , "\n")				
					}
		}
	#***
	# MAP
	maxval <- f.qk.yi[,1]
	indval <- 1
	for (tt in 2:TP){
		#tt <- 2
		m0 <- maxval
		maxval <- ifelse( f.qk.yi[,tt] > m0 , f.qk.yi[,tt] , m0 )
		indval <- ifelse( f.qk.yi[,tt] > m0 , tt , indval )	
				}
	MAP <- theta.k[ indval , ]
	#***
	# descriptives of classes
	
	score <- rowSums( dat2 *dat2.resp ) / rowSums( dat2.resp )
	plmat <- FALSE
	if ( is.vector(problevels) ){ 
		PL <- length(problevels) 
			}
	if ( is.matrix(problevels) ){ 
		PL <- nrow(problevels) 
		plmat <- TRUE
			}	
	theta.kk0 <- theta.k
	if ( model=="GOMRasch"){
	    PL <- 5
		problevels <- seq(0,1,len=PL)
		problevels2 <- c( problevels - diff(problevels)[1]/2 , 1.2 )
		for (kk in 1:K){
		# kk <- 1	
		for (ll in 1:PL){
	#		ll <- 1
			ind.kk <- which( ( theta.k[ , kk ] > problevels2[ll] ) & 
				( theta.k[ , kk ] <= problevels2[ll+1] ) )
			theta.kk0[ ind.kk , kk ] <- problevels[ll]
						} }
		pi.k <- pi.k[,1]					
			}  # end GOMRasch
		classdesc <- NULL
	if ( ! plmat ){				
		classdesc <- data.frame( matrix( 0 , 2*PL+1 , K ) )
		classdesc[1,] <- colSums( theta.kk0 * pi.k )
		for (kk in 1:K){
			for (ll in 1:PL){
				ll1 <- problevels[ll]
				classdesc[ ll+1 , kk ] <- sum( pi.k * ( theta.kk0[,kk] == ll1 ) )
				# aggregate( pi.k , list( theta.kk0[,kk]  ) , sum )[,2]
						}
				}
				
			rownames(classdesc)[1] <- "p.Class"
			colnames(classdesc) <- paste0("Class" , 1:K )
			rownames(classdesc)[2:(PL+1)] <- paste0( "p.problevel" , round( problevels , 3 ) , ".class" )
			rownames(classdesc)[PL+2:(PL+1)] <- paste0( "M.problevel" , round( problevels , 3 ) , ".class" )
		for (kk in 1:K){
			# kk <- 1
			ll.kk <- problevels
			for (ll in 1:PL){
				# ll <- 1
				ind.ll <- which( theta.kk0[,kk]  == problevels[ll] )
				rll.wt <- rowSums( f.qk.yi[ , ind.ll , drop=FALSE] )
				ll.kk[ll] <- weighted.mean( score , rll.wt )
					}
			classdesc[ PL + 2:(PL+1) , kk ]  <- ll.kk
				}
			}
	#***
	# output
	s2 <- Sys.time()
    res <- list("deviance" = dev , "ic"=ic , "item"=item , 
		"person" = person , "EAP.rel"=EAP.rel , 
		"MAP"=MAP , "classdesc" = classdesc , 
		"lambda"=lambda , "se.lambda"=se.lambda , 
		"mu"=mu , "Sigma"=Sigma , 
		"b"= b , "se.b" = se.b , 
		"f.yi.qk"=f.yi.qk , "f.qk.yi"=f.qk.yi , "probs"=probs ,
		"n.ik"=n.ik ,  "iter"=iter , "dat"=dat0 , "dat2"=dat2 , "dat2.resp"=dat2.resp , 
		"I"=I ,  "K"=K ,  "TP"=TP , "G"= 1 ,  "theta.k" = theta.k , "pi.k"=pi.k ,
		"problevels"=problevels , 	
		"model"=model , "s1"=s1 , "s2"=s2  , "plmat"= plmat
				)
	class(res) <- "gom"
	return(res)
		}
############################################################