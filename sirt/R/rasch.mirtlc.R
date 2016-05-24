
##########################################################
# Muiltidimensional Latent Class IRT models
rasch.mirtlc <- function( dat , Nclasses=NULL , modeltype="LC" , 
	dimensions = NULL , group = NULL ,
	weights=rep(1,nrow(dat)) , 
	theta.k = NULL , ref.item = NULL , 
	distribution.trait= FALSE ,  range.b =c(-8,8) , range.a =c(.2 , 6 ) ,  
	progress=TRUE , glob.conv=10^(-5)  ,
	conv1=10^(-5)  , mmliter=1000 , mstep.maxit = 3 , seed=0  , nstarts = 1 ,
	fac.iter = .35 ){
	#..................................................
	# preliminaries
	dat <- as.matrix(dat)
	theta.normal <- FALSE
	a <- c(1)
	# handle warnings
	warn_temp <- options()$warn
#	ref.item <- NULL
	if ( is.null(theta.k) ){ 
		  theta.fixed<-FALSE  } else { 
		  theta.fixed <- TRUE 
		  if( is.vector(theta.k) ){ Nclasses <- length(theta.k) }
		  if( is.matrix(theta.k) ){ Nclasses <- nrow(theta.k) }		  
			}
	s1 <- Sys.time()
	if ( is.null(group)){ group <- rep(1 , nrow(dat)) }
	G <- length( unique( group ))
	#	rescale weights
	
#    weights <- nrow(dat) * weights / sum(weights )
	dat2 <- dat
	dat2.resp <- 1* ( 1- is.na(dat ))
	dat2[ is.na(dat) ] <- 0
	# dat1: calculate matrix of frequencies 
	dat1 <- matrix( NA , nrow(dat) , Nclasses + 1 )
	dat1[,2] <- weights
	I <- ncol(dat2)     # number of items
	n <- nrow(dat2)     # number of persons	
	# initialize probabilities
	#******************
	if ( modeltype=="LC"){
		D <- 1
		pi.k <- rep( 1/Nclasses , Nclasses )
		if (G>1){ pi.k <- matrix( pi.k , nrow=Nclasses , ncol=G ) }
		lc.probs <- matrix( NA , nrow=Nclasses , ncol= I )
		b <- - stats::qnorm( colMeans( dat , na.rm=T ) )
		theta.k <- stats::qnorm( ( seq( 1 , Nclasses , 1 ) - .5 ) / Nclasses )
		for (cc in 1:Nclasses ){
			lc.probs[ cc , ] <- stats::plogis( theta.k[cc] - b )
					}
		pjk <- lc.probs
						}  # end LC
	#******************
	# MLC
	if ( modeltype %in% c("MLC1","MLC2") ){
		if ( is.null( dimensions ) ){ D <- 1 } else {
			D <- length( unique( dimensions ) )
			dimensions <- match( dimensions , unique( dimensions ) )
			distribution.trait <- FALSE
										}
		inut <- is.null(theta.k) 
	    if (inut){
			theta.k <- 2* stats::qnorm( seq( 1 / Nclasses / 2, 1 , 1/Nclasses ) )
							} else {
					if (D==1){ Nclasses <- length(theta.k )  }	# works for D=1
					if (D>1){ Nclasses <- nrow(theta.k )  }
							}

		if ( D==1){ pi.k <- dnorm( theta.k ) ; Qmatrix <- NULL }
		if ( D > 1 ){
			Qmatrix <- matrix( 0 , I , D )
			Qmatrix[ cbind( 1:I , dimensions ) ] <- 1
			if ( inut ){
				 theta.k <- 
						matrix( theta.k  , nrow=Nclasses , ncol=D , byrow=FALSE)
				 theta.kstart <- theta.k
				 if ( seed[1] != 0 ){
				 if (seed[1] >0){  set.seed( seed[1] ) } else{ set.seed( Sys.time() ) }
					 theta.k <- theta.k + matrix( rnorm( Nclasses*D , sd =2 ) , ncol=D )
								}
						} # end inut (is.null(theta.k))
            if ( D > 1 ){ Nclasses <- nrow( theta.k ) }
			pi.k <- rep(1/Nclasses , Nclasses )
								} # end if D > 1
		pi.k <- pi.k / sum( pi.k )
		b <- - stats::qlogis( colMeans( dat , na.rm=T ) )	
        a <- rep(1,I)
		if (G>1){ pi.k <- matrix( pi.k , nrow=Nclasses , ncol=G ) }
		lc.probs <- matrix( NA , nrow=Nclasses , ncol= I )
		if (D==1){
		for (cc in 1:Nclasses ){
			lc.probs[ cc , ] <- stats::plogis( theta.k[cc] - b )
					}
				}
		pjk <- lc.probs
		# design matrix theta
		des1 <- matrix( 0 , nrow=I*Nclasses , ncol=Nclasses )
		g1 <- cbind( seq(1,I*Nclasses) , rep( 1:Nclasses , each=I ) )
		des1[g1] <- 1
		des1 <- rbind( des1 , des1 )
		if ( D==1 ){ des.theta <- des1 }
		if ( D > 1 ){
			des.theta <- NULL
			for (dd in 1:D){
				# dd <- 1
				ind.dd <- which( dimensions != dd )
				sel.dd <- ( 1:( nrow(des1) ) ) %% I
				sel.dd[ sel.dd == 0 ] <- I
				des1.dd <- des1
				des1.dd[ sel.dd %in% ind.dd , ] <- 0
				des.theta <- cbind( des.theta , des1.dd )
							}
				}	
		# design matrix item difficulties
		des2 <- matrix( 0 , nrow=I*Nclasses , ncol=I )
		g1 <- cbind( seq(1,I*Nclasses) , rep( 1:I , Nclasses ) )
		des2[g1] <- -1
		des2 <- rbind( des2 , des2 )
		if ( theta.normal ){ theta.fixed <- TRUE }
		# select reference item		
		pval <- colMeans( dat , na.rm=T )
		p1 <- sort( pval , index.return=TRUE)$ix
		if ( D==1 & is.null(ref.item) ){ 
				ref.item <- p1[round(I/2)] 
						}
		if (D>1 & is.null(ref.item) ){ 
			ref.item <- NULL
			for (dd in 1:D){			
				# dd <- 1
				h1 <- which( dimensions == dd )
				ind.dd <- h1[ which.min( pval[h1] - 1/2 )[1] ]
				ref.item <- c( ref.item , ind.dd )
					}
					}
		
#		if ( theta.normal ){ ref.item <- NULL }
#		if ( ! is.null(ref.item ) ){ 	
		des2 <- des2[ , - ref.item ]#		}
		# reference item
		b[ ref.item ] <- 0
		est.a <- seq(1,I)
		est.a[ref.item] <- 0
#		des.theta <- des1
		des.b <- des2
						}  # end MLC	
	#***************
		# iterate over different starts
		if ( ( seed[1] < 0) | ( length(seed) < nstarts) ){ 
			seed <- round( stats::runif( nstarts , 1 , 10000 )) 
					}
		devL <- rep(NA , nstarts )
		NN1dev <-  1*10^90
		#####################################
		# different starts
		for (nn in 1:nstarts ){ 
			seed.nn <- seed[nn]
			if ( seed[nn] > 0 ){
				set.seed( seed.nn )
				pjk <- matrix( stats::runif( I*Nclasses ) , nrow= Nclasses , ncol=I )
#				theta.k <- theta.k + matrix( rnorm( Nclasses*D , sd =.7 ) , ncol=D 
#				theta.k <- theta.kstart + matrix( rnorm( Nclasses*D , sd =.07 ) , ncol=D )
			if (modeltype%in%c("MLC1","MLC2") & ( nstarts > 1 ) & ( ! theta.fixed ) & (D>1) ){ 
					pi.k <- stats::runif( Nclasses , 0 , 1 ) 
					pi.k <- pi.k / sum( pi.k )
					theta.k <- .7*theta.k + matrix( stats::rnorm( Nclasses*D , sd =.97 ) , ncol=D )
								}
			if (modeltype%in%c("MLC1","MLC2") & ( nstarts > 1 ) & ( ! theta.fixed ) & (D==1) ){ 
					pi.k <- stats::runif( Nclasses , 0 , 1 ) 
					pi.k <- pi.k / sum( pi.k )
					theta.k <- .7*theta.k + stats::rnorm( Nclasses*D , sd =.97 ) 
								}
								}
	dev <- 1
	iter <- 0
	dev.change <- par.change <- 1000
	# display
	disp <- "...........................................................\n"   
    NN1dev <-  1*10^90
	
	
	######################################
	# begin EM algorithm
	while ( ( dev.change > glob.conv | par.change > conv1 ) & iter < mmliter ){
			if (progress){ 
			  cat(disp) 
			  cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )
			  flush.console()
						}   
	pjk0 <- pjk
	pi.k0 <- pi.k
	theta.k0 <- theta.k
	dev0 <- dev
# 	b0 <- b
# a0 <- Sys.time()
	# E step latent class analysis
	if ( modeltype =="LC"){
		res1 <- .e.step.mirtlc.lc( dat1 , dat2 , dat2.resp , pi.k , pjk , I , 
						   group , G ,  theta.k , f.qk.yi=NULL  )
							}
	if ( modeltype %in% c("MLC1","MLC2")){
		res1 <- .e.step.mirtlc.mlc1( dat1 , dat2 , dat2.resp , pi.k , pjk , I , 
						   b , a , group , G ,  theta.k , 
						   D , dimensions , Qmatrix , f.qk.yi=NULL  )
							}
# a1 <- Sys.time() ; adiff <- a1-a0 ; cat("\ne step" , adiff ) ; a0 <- a1							
						
	n.k <- res1$n.k ; n.jk <- res1$n.jk 
	r.jk <-  res1$r.jk ; f.qk.yi <- res1$f.qk.yi
	pjk <- res1$pjk   ; f.yi.qk <- res1$f.yi.qk
	ll <- res1$ll   ; dev <- -2*ll
	# Mstep:
	if ( modeltype =="LC"){
		res2 <- .m.step.mirtlc.lc( pjk , n.k , r.jk , n.jk , G , Nclasses)
					      }
	if ( modeltype%in%c("MLC1","MLC2") ){
	options(warn=-1) 
		res2 <- .m.step.mirtlc.mlc1( pjk , n.k , r.jk , n.jk , G , Nclasses ,
				theta.k , b , a , I , ref.item , mstep.maxit ,
				des.theta , des.b , theta.fixed , theta.normal ,
				f.qk.yi	, D  , distribution.trait , est.a , Qmatrix , modeltype,
				 range.b , range.a , iter , fac.iter )
		b <- res2$b ; theta.k <- res2$theta.k ; pi.k <- res2$pi.k
		a <- res2$a
	options(warn=warn_temp)	
					        }					  
							
													
# a1 <- Sys.time() ; adiff <- a1-a0 ; cat("\nm step" , adiff ) ; a0 <- a1
	pi.k <- res2$pi.k ;  pjk <- res2$pjk
	# prevent label switching (modeltype == "LC")
	if ( modeltype == "LC"){
		ind <- order( rowMeans( pjk ) )
		pjk <- pjk[ ind , , drop=FALSE]
		if (G==1){ pi.k <- pi.k[ ind ] }
		if (G>1){ pi.k <- pi.k[ ind , ] }
							}

							
							
	# convergence criteria
	a1 <- max( abs( pjk - pjk0 ) )
	a2 <- max( abs( pi.k - pi.k0 ) )
	if ( modeltype%in%c("MLC1","MLC2") ){ 
			a3 <- max( abs( theta.k - theta.k0 ))  } else { a3 <- 0 }
	dev.change <- abs( ( dev - dev0)/ dev0 )
	par.change <- max( c(a1,a2,a3))
	iter <- iter + 1
	
	# settings
	if ( dev < NN1dev ){
		NN1pjk <- pjk ; 		NN1pi.k <- pi.k
		NN1dev <- dev ; 		NN1ll <- ll
		NN1res1 <- res1 ; NN1theta.k <- theta.k
		NN1a <- a ; NN1b <- b ; NN1iter <- iter			
						}
						
		if ( progress  ){   
				   cat( paste( "   Deviance = "  , 
					round( dev , 4 ) , 
					if (iter > 0 ){ " | Deviance change = " } else {""} ,
					if( iter>0){round( - dev + dev0 , 6 )} else { ""}  ,
					 "  START " , nn , " (Seed " , seed.nn , ")\n",sep="") )
					# maximum probability change       
					cat( paste( "    Maximum item parameter change = " , 
							paste( round(a1 ,6) , collapse=" " ) , "\n" , sep=""))
					# maximum probability distribution parameter change
					cat( paste( "    Maximum probability distribution change = " , 
							paste( round(a2 ,6) , collapse=" " ) , "\n" , sep=""))
					cat( paste( "    Maximum theta parameter change = " , 
							paste( round(a3 ,6) , collapse=" " ) , "\n" , sep=""))
							}               
	#    print( paste( iter , ll) ) ; flush.console()
		}
	##### end algorithm  #***************************	

		NN1pjk -> pjk ; 		NN1pi.k -> pi.k
		NN1dev -> dev ; 		NN1ll -> ll
		NN1res1 -> res1 ; NN1theta.k -> theta.k
		NN1a -> a ; NN1b -> b ; NN1iter -> iter			
	
	##############################################################	
	# collect results of nstarts
	if ( ( nn == 1 ) ){
		NNpjk <- pjk ; 		NNpi.k <- pi.k
		NNdev <- dev ; 		NNll <- ll ;
		NNres1 <- res1 ; NNtheta.k <- theta.k
		NNa <- a ; NNb <- b ; NNiter <- iter
					}
	devL[nn] <- dev
	if ( ( nn > 1 ) & ( dev < NNdev ) ){
		NNpjk <- pjk ; 		NNpi.k <- pi.k
		NNdev <- dev ; 		NNll <- ll
		NNres1 <- res1 ; NNtheta.k <- theta.k
		NNa <- a ; NNb <- b ; NNiter <- iter
					}
					
				}
				
	if (nstarts > 1){ 
			NNpjk -> pjk ; 		NNpi.k -> pi.k
			NNdev -> dev ; 		NNll -> ll
			NNres1 -> res1 ; NNtheta.k -> theta.k
			NNa -> a ; NNb -> b ; NNiter -> iter
					}
	#############################################
	# labels for pjk
	colnames(pjk) <- colnames(dat)
	rownames(pjk) <- paste("Class" , 1:Nclasses , sep="")
	if (G==1 ){ names(pi.k) <- rownames(pjk) }
	if (G>1 ){ rownames(pi.k) <- rownames(pjk) }	
		
	# Information criteria
    # calculations for information criteria
    ic <- list( "deviance" = dev , "n" = nrow(dat) )
    # number of parameters to be estimated
	if ( modeltype == "LC" ){
		ic[[ "itempars" ]] <- Nclasses * I
		ic[[ "traitpars" ]] <- ( Nclasses - 1 ) * G
						}
	if ( modeltype%in%c("MLC1","MLC2") ){
		ic[[ "itempars" ]] <- I - D
		if ( modeltype=="MLC2"){ ic$itempars <- 2*(I-D) }
		ic[[ "traitpars" ]] <- D*Nclasses + ( Nclasses - 1 ) * G
		# trait distribution + probabilities
		ic$traitpars <- ic$traitpars - D*theta.fixed*Nclasses
		if ( distribution.trait == "normal" ){
			ic$traitpars <- 2*G
										}		
		if ( distribution.trait == "smooth2" ){	ic$traitpars <- 2*G	}		
		if ( distribution.trait == "smooth3" ){	ic$traitpars <- 3*G	}		
		if ( distribution.trait == "smooth4" ){	ic$traitpars <- 4*G	}		
										}

	ic$np <- ic$itempars + ic$traitpars
	
	ic$n <- n # number of persons
	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC (conistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	
	###########################
	# item and trait parameter
	cor.trait <- skewness.trait <- mean.trait <- 
			sd.trait <- trait <- item <- NULL
	if (modeltype%in%c("MLC1","MLC2") ){
		mean.trait <- rep(NA,G)
		if (D>1){ skewness.trait <- sd.trait <- mean.trait <- matrix(NA , D , G ) }
		item <- data.frame("b" = b , "a" = a )
		rownames(item) <- colnames(dat)
		item$thresh <- b*a
		
		trait <- data.frame( "theta" = theta.k )
		if ( D>1){ colnames(trait) <- paste( "theta.Dim" , 1:D , sep="") }
		# mean
		if (D==1){ 
			mean.trait[1] <- M <- sum( pi.k[,1] * theta.k )
#			item$b.cent <- item$b - M
#			trait$theta.cent <- theta.k - M
			mean.trait <- sapply(1:G , FUN = function(gg){ 
									sum( pi.k[,gg] * theta.k )  } )
			sd.trait <- sqrt( sapply(1:G , FUN = function(gg){ 
							sum( pi.k[,gg] * theta.k^2 )  } ) - mean.trait^2 )

			skewness.trait <- sapply(1:G , FUN = function(gg){ 
							m3 <- sum( pi.k[,gg] * ( theta.k - mean.trait[gg] )^3 )
							m3 / sd.trait[gg]^3
								} )
						}
			# mean and SD trait in multidimensional case
    		if (D> 1){
			    cor.trait <- array( NA , dim=list(D,D,G) )
				for (dd in 1:D){
				for (gg in 1:G){
#					dd <- 1 ; 				gg <- 1
					mean.trait[ dd , gg ] <- sum( theta.k[,dd] * pi.k[,gg] )
					sd.trait[dd,gg] <- sqrt( sum( theta.k[,dd]^2 * pi.k[,gg] ) - mean.trait[dd,gg]^2 )
     				m3 <- sum( pi.k[,gg] * ( theta.k[,dd] - mean.trait[gg] )^3 )
					skewness.trait[dd,gg] <- m3 / sd.trait[gg]^3
								} 
							}
				for (dd1 in 1:D){
				for (dd2 in 1:D){
				for (gg in 1:G){
#					dd <- 1 ; 				gg <- 1
					mdd1 <- mean.trait[ dd1 , gg ]
					mdd2 <- mean.trait[dd2,gg]
					sdd1 <- sd.trait[ dd1 , gg ]
					sdd2 <- sd.trait[dd2,gg]
					cor.trait[dd1,dd2,gg] <- sum( theta.k[,dd1] * theta.k[,dd2] * pi.k[,gg] )
					cor.trait[dd1,dd2,gg] <- ( cor.trait[dd1,dd2,gg] - mdd1 * mdd2 ) / ( sdd1 * sdd2 )
								}
							} }							
							
						} # end D>1
			if (G==1){ trait$prob <- pi.k }
			if (G>1){ trait <- cbind( trait , pi.k ) }
			trait <- data.frame( trait )
			if (G>1){ 
				colnames(trait)[ seq( ncol(trait) -G+1 , ncol(trait) , 1) ] <- 
							paste("prob.Group" , 1:G , sep="") }
			# standardized theta dimensions			
			if (D>1){
			 for (dd in 1:D){
#			dd <- 1
				trait[ , paste("theta.stand.Dim" , dd,sep="") ] <- 
						( trait[ , dd ] - mean.trait[dd,1] ) / sd.trait[dd,1]
							} # end dd in 1:D
						} # end D>1
					}  # end MLC1
	#####################################
    # output LC
	if ( modeltype == "LC"){	
		theta.k <- diag( length( theta.k) )
							}
							
	##################################################
	# item response probabilities
	d1 <- dim(pjk)
	rprobs <- array( 0 , dim=c( d1[2] , 2 , d1[1] ) )
	dimnames(rprobs)[[1]] <- colnames(dat)
	rprobs[,2,] <- t( pjk )
	rprobs[,1,] <- 1 - t(pjk)
							
					
	###########################
	# collect results
	s2 <- Sys.time()
	res <- list( "pjk" = pjk , "rprobs" = rprobs , "pi.k" = pi.k , "theta.k" = theta.k ,
			"item" = item , "trait" = trait , "mean.trait" = mean.trait ,
			"sd.trait" = sd.trait , "skewness.trait" = skewness.trait , 
			"cor.trait" = cor.trait , 
			"ic" = ic , "D" = D , "G" = G , 
			"deviance" = dev , "ll" = ll , "Nclasses" = Nclasses , 
			"modeltype"=modeltype , "estep.res" = res1 , "dat" = dat , 
			"devL" = devL , "seedL" = seed ,  "iter"=iter , "s1" = s1 , "s2" = s2 ,
			"distribution.trait"=distribution.trait)
	class(res) <- "rasch.mirtlc"
	return(res)		
		}
####################################################################






