

###########################################
# Structured latent class analysis
###########################################

slca <- function( data , group=NULL, 
    weights=rep(1, nrow(data)), 
	Xdes ,  Xlambda.init=NULL , Xlambda.fixed=NULL , 
	Xlambda.constr.V=NULL , Xlambda.constr.c=NULL , 
	delta.designmatrix =NULL ,  delta.init = NULL , 
	delta.fixed = NULL , delta.linkfct = "log" ,  
    maxiter=1000, conv=10^(-5), globconv=10^(-5), msteps=4 , 
	convM=.0005 , decrease.increments = FALSE , oldfac = 0 , 
	seed=NULL , progress = TRUE ,  ...){	
	#************************************************************
	# mean.constraint [ dimension , group , value ]
	# Sigma.constraint [ dimension1 , dimension2 , group , value ]	
	cl <- match.call()
    theta.k <- NULL	
	#*************************
	# data preparation
	s1 <- Sys.time()
	e1 <- environment()	
	use.freqpatt <- FALSE 
	## prevent from warnings in R CMD check "no visible binding"
	## gdm: no visible binding for global variable 'TD'
	TD <- TP <- EAP.rel <- mean.trait <- sd.trait <- skewness.trait <- NULL
	K.item <- correlation.trait <- NULL 
    se.theta.k <- NULL	
	data0 <- data <- as.matrix(data)
	dat.resp0 <- dat.resp <- 1 - is.na(data)
	dat <- data
	dat[ is.na(data) ] <- 0
	dat0 <- dat

	# maximal categories
	K <- max(dat)
	maxK <- K+1
	# list of indicator data frames
	dat.ind <- as.list( 1:(K+1) )
	for (ii in 0:K){
		dat.ind[[ii+1]] <- 1 * ( dat==ii )*dat.resp
				}
	I <- ncol(dat)	# number of items
	n <- nrow(dat)
# cat("a300\n")	
	# arrange groups
	if ( is.null(group)){ 
		G <- 1 
		group0 <- group <- rep(1,n)
				} 
			else {
		group0 <- group
	  if( ! ( is.numeric(group) ) ){
			gr2 <- unique( sort(paste( group ) ))
			    } else {
		gr2 <- unique( sort( group ) )
						}
		G <- length(gr2)
		group <- match( group , gr2 )
							}
	group.stat <- NULL
	if (G>1){
		a1 <- stats::aggregate( 1+0*group , list(group) , sum )
		a2 <- rep("",G)
		for (gg in 1:G){
			a2[gg] <- group0[ which( group == gg )[1]  ]
						}
		group.stat <- cbind( a2 , a1 )
		colnames(group.stat) <- c(  "group.orig" , "group" , "N"  )
	    Ngroup <- a1[,2]		
			}	
    if (G==1){ Ngroup <- length(group) }
# cat("a400\n")
	KK <- K	# if KK == 1 then a slope parameter for all items is estimated
    deltaNULL <- 0
	if ( is.null(delta.designmatrix) ){
	    deltaNULL <- 1
		delta.designmatrix <- diag( dim(Xdes)[3] )
							}
# print(delta.designmatrix)	
	# lambda basis parameters for X
	if ( is.null(seed) ){
		seed.used <- round( stats::runif(1,2,10000 ) )
					} else { seed.used <- seed }
	set.seed( seed.used )
	
	Nlam <- dim(Xdes)[[4]]
	if ( is.null( Xlambda.init ) ){
		Xlambda.init <- stats::runif( Nlam , -1 , 1 )
						}
	Xlambda <- Xlambda.init
	
	# starting values for distributions
	
#    delta.init <- matrix(1,nrow=1,ncol=1)
	TP <- nrow(delta.designmatrix)

	if (  ! is.null(delta.init) ){	
		delta <- delta.init
		if ( delta.linkfct == "log"){
			pik <- exp( delta.designmatrix %*% delta.init[,1] )
						} else {
			pik <- stats::plogis( delta.designmatrix %*% delta.init[,1] )
							}
									} else {
		# delta.init <- runif( Nlam , -.5 , .5 )
		pik <- rep( 1 /TP , TP ) + runif(TP , 0 , .5 )
		pik <- pik / sum(pik)
		if ( delta.linkfct == "logit"){
		g1 <- solve( crossprod( delta.designmatrix )) %*% 
					t( delta.designmatrix) %*% pik				
				delta <- matrix( g1[,1] , 
						nrow=ncol(delta.designmatrix) , ncol=G)
						} else {
			delta <- matrix( 0 , ncol(delta.designmatrix) , G )		
		delta[1,] <- 1		}
					}

	pi.k <- matrix( 0 , TP , G )
	for (gg in 1:G){ pi.k[,gg] <- pik }
	n.ik <- array( 0 , dim=c(TP,I,K+1,G) )	


	#***
	# response patterns
    resp.ind.list <- list( 1:I )
	for (ii in 1:I){ resp.ind.list[[ii]] <- which( dat.resp[,ii] == 1)  }	
#cat("a900\n")	

	se.Xlambda <- 0*Xlambda
	max.increment.Xlambda <- .3

	# lambda constraints
	Xlambda.constraint <- NULL


	#***
	# preparations for calc.counts
	dat.ind2 <- as.list( 1:(K+1) )
	ind.group <- as.list( 1:G )
    for (kk in 1:(K+1)){ 
		l1 <- as.list(1:G)
		for (gg in 1:G){
		  if ( ! use.freqpatt ){
			ind.gg <- which( group == gg )
			ind.group[[gg]] <- ind.gg
			dkk <- (dat.ind[[kk]])[ ind.gg , ]
			l1[[gg]] <- dkk * dat.resp[ind.gg,] * weights[ind.gg] 	
							}
		  if ( use.freqpatt ){
			dkk <- dat.ind[[kk]]
			if (G>1){ 	wgg <- weights[,gg]	 }
			if (G==1){ wgg <- weights 
					ind.group[[gg]] <- which( group==gg)
							}
			l1[[gg]] <- dkk * dat.resp * wgg 
							}
						}   # end gg
			dat.ind2[[kk]] <- l1
					}

# print("p400")					
	
	#*****
	# reducing computational burden for design matrix
	dimXdes <- dim(Xdes)
	Xdes_ <- as.vector(Xdes)
#	res <- calc_Xdes( Xdes_ , dimXdes )			# Rcpp function  .Call
	res <- .Call( "calc_Xdes" , Xdes_ , dimXdes , package="CDM" )
	# XdesM     [ii,kk,tt,ll, value ] 
	NX <- res$NXdesM
	XdesM <- res$XdesM[1:NX,]
	XdesM <- XdesM[ order( XdesM[,1]*NX + XdesM[,3] ) , ]	
	
	# Xlambda constraints	
	if ( ! is.null(Xlambda.constr.V) ){
			V <- Xlambda.constr.V
			e2 <- matrix( Xlambda.constr.c , nrow=ncol(V) , ncol=1 )
			V1 <- solve( crossprod(V) )
										}				
					
	#---
	# initial values algorithm
	dev <- 0	; iter <- 0
	globconv1 <- conv1 <- 1000
	disp <- paste( paste( rep(".", 70 ) , collapse="") ,"\n", sep="")
	mindev <- Inf
# cat("s700\n")	
	############################################
	# BEGIN MML Algorithm
	############################################
		
	while( ( iter < maxiter ) & ( ( globconv1 > globconv) | ( conv1 > conv) ) ){
		
		#****
		# collect old parameters
		Xlambda0 <- Xlambda 
		dev0 <- dev
		delta0 <- delta
		pi.k0 <- pi.k

#  z0 <- Sys.time()
 		
		#****
		#1 calculate probabilities
#		probs <- .slca.calc.prob( Xdes , Xlambda , I,K,TP)
		probs <- .slca.calc.prob( XdesM , dimXdes , Xlambda )

#  cat("calc.prob") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		


		#*****
		#2 calculate individual likelihood
		h1 <- matrix( 1 , nrow=n , ncol=TP )
		res.hwt <- calc_posterior.v2(rprobs= probs , gwt=h1 , 
					 resp=dat , nitems= I , 
					 resp.ind.list=resp.ind.list , normalization=FALSE , 
					 thetasamp.density= NULL , snodes=0 )	
		p.xi.aj <- res.hwt$hwt 	

# cat("calc.like") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
		
		#*****
		#3 calculate posterior and marginal distributions
		res <- .gdm.calc.post(pi.k,group,p.xi.aj,weights,G,ind.group,
				use.freqpatt )
		p.aj.xi <- res$p.aj.xi
		pi.k <- res$pi.k


# cat("calc.post") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
		
		#*****
		#4 calculate expected counts
		# n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
		res <- .slca.calc.counts(G,weights,dat.ind,dat,
					dat.resp,p.aj.xi,K,n.ik,TP,I,group , dat.ind2 , ind.group ,
					use.freqpatt )		
		n.ik <- res$n.ik
		n.ik1 <- res$n.ik1
		N.ik <- res$N.ik
		N.ik1 <- res$N.ik1

# cat("calc.counts") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		


		#*****
		#5 M step: Xdelta parameter estimation
		# n.ik  [1:TP,1:I,1:K,1:G]
		# probs[1:I,1:K,1:TP]
		res <- .slca.est.Xlambda(Xlambda , Xdes , probs, n.ik1, N.ik1, I, K, G,
			max.increment,TP,msteps,convM , Xlambda.fixed , XdesM , dimXdes , oldfac)		
	
		Xlambda <- res$Xlambda
		se.Xlambda <- res$se.Xlambda
		if (decrease.increments){ 	
				max.increment.Xlambda <- res$max.increment.Xlambda/1.01	
				max.increment <- max.increment.Xlambda 
							}
		
		# linear constraints on Xlambda parameters
		
		# below is code copied from rasch.pml3 (sirt package)
		#................
			# linear constraints: Let e be the vector of error
			# correlations, V a design matrix and c a vector.
			# The constraints can be written in the form 
			# c = V * e . Then V*e - c = 0.
			# See the Neuhaus paper:
			# e_cons = e + V * (V'V)^(-1) * ( c - V * e )
		if ( ! is.null(Xlambda.constr.V) ){
            e1 <- matrix( Xlambda , ncol=1 )
			Xlambda <- ( e1 + V %*% V1 %*% ( e2 - t(V) %*% e1 ) )[,1]		
									}
		#	if ( ! is.null( err.constraintM ) ){
		#	    V1 <- solve( t(V) %*% V )
		#		e1 <- matrix( eps.corr , ncol=1 )
		#		eps.corr2 <- eps.corr + V %*% V1 %*% ( err.constraintV - t(V) %*% e1 ) 
		#		eps.corr <- eps.corr2
		#								}		
		#    err.constraintM=NULL , err.constraintV=NULL ,
		#     V <- err.constraintM
		#     err.constraintV <- ...
		#     Me=v
		#.................
		#          
		#     c = V*Xlambda
		
		
		# Xlambda.constr.V=NULL , Xlambda.constr.c=NULL , 



# cat("est b") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		
	
						
		#*****
		#7 M step: estimate reduced skillspace
		res <- .slca.est.skillspace(Ngroup, pi.k , 
			delta.designmatrix , G , delta  , delta.fixed , 
			eps=10^(-4) , oldfac , delta.linkfct)
		pi.k <- res$pi.k
		delta <- res$delta
		covdelta <- res$covdelta

# cat("skillspace") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		

			
		#*****
		#8 calculate likelihood
		# n.ik [ TP , I , K+1 , G ]
		# N.ik [ TP , I , G ]
		# probs [I , K+1 , TP ]
		ll <- 0
		for (gg in 1:G){ 
#			ind.gg <- which(group==gg)
			if ( ! use.freqpatt ){
				ind.gg <- ind.group[[gg]]
				ll <- ll + sum( weights[ind.gg] * log( rowSums( p.xi.aj[ind.gg,] * 
							matrix( pi.k[,gg] , nrow= length(ind.gg) , ncol=nrow(pi.k) , byrow=TRUE ) ) ) )
									}
			if ( use.freqpatt ){
			  if (G>1){   wgg <- weights[,gg]  }
			  if (G==1){ wgg <- weights }
				ll <- ll + sum( wgg * log( rowSums( p.xi.aj * 
							matrix( pi.k[,gg] , nrow= nrow(p.xi.aj) , 
									ncol=nrow(pi.k) , byrow=TRUE ) ) ) )
									}
							
							}
		dev <- -2*ll	

#  cat("calc LL") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		

		#~~~~~~~~~~~~~~~~~~~~~~~~~~
		# display progress		
		gg0 <- 0
		gg1 <- abs( Xlambda - Xlambda0 )
# print(Xlambda)		
		pardiff <- max( gg1 , gg0 )
		deltadiff <- abs( pi.k - pi.k0 )	
		conv1 <- max( c(pardiff,deltadiff))
		globconv1 <- abs( dev - dev0) 
		iter <- iter +1
		if (progress){
			cat(disp)	
			cat("Iteration" , iter , "   " , paste( Sys.time() ) , "\n" )		
			cat( paste( "   Deviance = "  , 
					   round( dev , 4 ) , 
					  if (iter > 1 ){ " | Deviance change = " } else {""} ,
							 if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,sep=""))
			if ( dev > dev0 & (iter>1 ) ){ cat( "  Deviance increases!") } ; cat("\n")
			cat( paste( "    Maximum Xlambda parameter change = " , 
								 round( max( gg1) , 6 ) ,  " \n"   )  )  
			cat( paste( "    Maximum distribution parameter change = " , 
								 round( max( deltadiff ) , 6 ) ,  " \n"   )  )  
			utils::flush.console()
				}

		# save values corresponding to minimal deviance
		if ( ( dev < mindev ) | ( iter == 1 ) ){
			Xlambda.min <- Xlambda
			se.Xlambda.min <- se.Xlambda
			pi.k.min <- pi.k
			n.ik.min <- n.ik
			probs.min <- probs
			delta.min <- delta
			covdelta.min <- covdelta
			mindev <- dev
			iter.min <- iter
					}
		
								}
		############################################
		# END MML Algorithm
		############################################
		
			Xlambda.min -> Xlambda
			se.Xlambda.min -> se.Xlambda
			pi.k.min -> pi.k
			n.ik.min -> n.ik
			probs.min -> probs
			delta.min -> delta
			covdelta.min -> covdelta
			mindev -> dev
			# iter.min -> iter
		
		# names
		if ( is.null(dimnames(Xdes)[[4]] ) ){
			dimnames(Xdes)[[4]] <- paste0("lam" , 1:Nlam ) }
		if ( is.null(dimnames(Xdes)[[3]] ) ){
			dimnames(Xdes)[[3]] <- paste0("Class" , 1:TP ) }	
		
		names(Xlambda) <- dimnames(Xdes)[[4]]
		colnames(pi.k) <- paste0("Group" , 1:G )
		rownames(pi.k) <- dimnames(Xdes)[[3]]		
		
		#  collect item parameters
		item1 <- array( aperm( probs , c(2,1,3)) , dim= c(I*maxK , TP) )
		colnames(item1) <- dimnames(Xdes)[[3]] 
		item <- data.frame("item" = rep(colnames(dat) , each=maxK) , 
						"Cat" = rep(0:K , I) , item1 )			
		rownames(item) <- paste0( rep(colnames(dat) , each=maxK) , "_Cat" , rep(0:K , I) )
		
		
#		res <- .gdm.collect.itempars( data , K , D , b , a , TD , thetaDes ,
#					irtmodel , se.b , se.a , data0)
#		.attach.environment( res , envir=e1 )					
				
		# calculate distribution properties
#		res <- .gdm.calc.distributionmoments( D , G , pi.k , theta.k )
#		.attach.environment( res , envir=e1 )	
		
		# Information criteria
		ic <- .slca.calc.ic( dev , dat , G ,   
					K, TP ,I , delta.designmatrix , delta.fixed ,
					Xlambda , Xlambda.fixed , data0 , deltaNULL ,
					Xlambda.constr.V
						)

	#########################################
	# item fit [ items , theta , categories ] 
	# # n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	probs <- aperm( probs , c(3,1,2) )
#	itemfit.rmsea <- itemfit.rmsea( n.ik , pi.k , probs ,
#			itemnames = colnames(data) )		
#   item$itemfit.rmsea <- itemfit.rmsea$rmsea

	# person parameters
#	res <- .gdm.person.parameters( data , D , theta.k , p.xi.aj , p.aj.xi , weights )	
#	.attach.environment( res , envir=e1 )
	
	# person parameters
	mle.class <- max.col( p.xi.aj ) 
	map.class <- max.col( p.aj.xi ) 
	
	
	
	#*************************
	# collect output	
	s2 <- Sys.time()
	res <- list( "item" = item , # "person" = person , 
				"deviance"=dev , "ic" = ic , 
				"Xlambda" = Xlambda , "se.Xlambda" = se.Xlambda , 
				"pi.k"=pi.k , 	
				"pjk" = probs , "n.ik" = n.ik ,  
				"G"=G , "I" = ncol(data) , "N" = nrow(data) , 
				"TP"=TP , 
				"delta" = delta , "covdelta"=covdelta , 
				"delta.designmatrix" = delta.designmatrix , 
				"MLE.class"=mle.class , "MAP.class" = map.class , 
				"data" = data ,
				"group.stat"=group.stat )
	res$p.xi.aj <- p.xi.aj ; res$posterior <- p.aj.xi 
	res$K.item <- K.item
	res$time <- list( "s1"=s1,"s2"=s2 , "timediff"=s2-s1)
	res$iter <- iter
	res$iter.min <- iter.min
	res$converged <- iter < maxiter
	
	# some further values for modelfit.gdm
	res$AIC <- res$ic$AIC
	res$BIC <- res$ic$BIC	
	res$Npars <- res$ic$np	
	res$loglike <- - res$deviance / 2
	res$seed.used <- seed.used
	res$Xlambda.init <- Xlambda.init
	res$delta.init <- delta.init
	
	res$control$weights <- weights
	res$control$group <- group
	
	if (progress){
                cat("----------------------------------- \n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("----------------------------------- \n")
				}
	res$call <- cl			
	class(res) <- "slca"
	return(res)
				
		}		
###################################################
