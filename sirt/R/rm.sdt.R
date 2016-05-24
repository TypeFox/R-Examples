
#################################################################
# Hierarchical rater model
# MML estimation
rm.sdt <- function( dat , pid , rater ,Qmatrix=NULL , theta.k=seq(-9,9,len=30) , 	
	est.a.item=FALSE , est.c.rater= "n" , 
	est.d.rater= "n" , est.mean = FALSE , skillspace="normal" , 
	tau.item.fixed = NULL , a.item.fixed = NULL , 
	d.min=.5 , d.max=100 ,  d.start = 3 , 
	max.increment=1 , numdiff.parm=.00001 , maxdevchange=.10 ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001){
	#..........................................................
	s1 <- Sys.time()
	theta.k0 <- theta.k
	max.b.increment <- max.increment
	pi.k <- dnorm( theta.k )
	pi.k <- pi.k / sum( pi.k )
	# process data
	procdata <- res <- .prep.data.rm( dat=dat , rater=rater , pid=pid )
	dat2 <- as.matrix(res$dat2)
	dat2.resp <- as.matrix(res$dat2.resp)
	rater.index1 <- res$rater.index
	dataproc.vars <- res$dataproc.vars
	VV <- res$VV
	RR <- res$RR
	item.index <- res$dataproc.vars$item.index 
	rater.index <- res$dataproc.vars$rater.index 
	
	# maximum categories
	maxK <- apply( dat , 2 , max , na.rm=T )
	K <- max( maxK )
	if ( is.null(Qmatrix) ){
		Qmatrix <- matrix( 1:K , nrow=VV , ncol=K , byrow=TRUE)
					}
	TP <- length(theta.k)
	I <- VV*RR
	
	# define constraints on tau.item parameters
	# if not all categories are observed
#	tau.item.fixed <- NULL
	if ( is.null( tau.item.fixed )){
		if ( min(maxK) < K ){
			for (vv in 1:VV){
			#vv <- 1
				K.vv <- maxK[vv]
				if ( K.vv < K ){
					for (zz in (K.vv+1):K ){
						d1 <- data.frame( "item"= vv, "categ"=zz , "val"=99 )	
						tau.item.fixed <- rbind( tau.item.fixed , d1 ) 				
										}
								}
						}
			tau.item.fixed <- as.matrix(tau.item.fixed )
					}
				}
				
	# starting values for item difficulties
	tau.item <- matrix( 0 , nrow=VV , ncol=K )
	rownames(tau.item) <- colnames(dat)
	
	tau.item <- matrix( seq( -2 , 2 , len=K ) , nrow=VV , ncol=K , byrow=TRUE )
	if ( ! is.null(tau.item.fixed) ){
	    tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
									}
		
	a.item <- rep( 1 , VV )
	if (skillspace == "discrete" ){
			est.mean <- TRUE
						}
	if ( ! is.null( a.item.fixed ) ){
		est.a.item <- TRUE
		a.item[ a.item.fixed[,1] ] <- a.item.fixed[,2]
									}
	
	# rater parameter
	d.rater <- matrix( d.start , nrow=I , ncol=1 )
	c.rater <- matrix( d.start*((1:K) - .5 ) , nrow=I , ncol=K , byrow=TRUE )
	
	# set c.rater for fixed items to 99
	c.rater.fixed <- NULL
	if ( ! is.null( tau.item.fixed ) ){
		tau1 <- tau.item.fixed[ tau.item.fixed[,3] == 99 , , drop=FALSE]		
		ind <- match( item.index , tau1[,1] )
		c.rater.fixed <- tau1[ ind , ]
        c.rater.fixed[,1] <- seq( 1 , nrow(c.rater.fixed) )	
#		c.rater.fixed[,3] <- - 99
		c.rater.fixed[,3] <- 999
		c.rater.fixed <- c.rater.fixed[ ! is.na( c.rater.fixed[,2] ) , ] 
		c.rater[ c.rater.fixed[,1:2] ] <- c.rater.fixed[,3]
		
							}
		
    # init standard errors
	se.d.rater <- NA*d.rater
	se.c.rater <- NA*c.rater
	se.a.item <- NA*a.item
		
	# inits
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	sigma <- 1
	disp <- "...........................................................\n"	
	prob.item <- NULL ; prob.rater <- NULL
	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter )
						){
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
		
		# previous values
		d.rater0 <- d.rater
		tau.item0 <- tau.item
		dev0 <- dev
		sigma0 <- sigma
		a.item0 <- a.item
		c.rater0 <- c.rater

  a0 <- Sys.time()		
		# calculate probabilities
		res <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				theta.k ,RR )		
# cat("calcprob "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;			

        probs <- res$prob.total				
		prob.rater <- res$prob.rater
		prob.item <- res$prob.item
		
		# calculate posterior
		res <- .rm.posterior( dat2 , dat2.resp , TP , pi.k , K, I , probs )
		f.yi.qk <- res$f.yi.qk
		f.qk.yi <- res$f.qk.yi
		n.ik <- res$n.ik
		N.ik <- res$N.ik
		pi.k <- res$pi.k
		ll <- res$ll

# cat("posterior "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;			
																					
    	# estimate tau.item parameters
		if (iter ==0){	max.b.increment -> tau.item.incr }
		res <- .rm.hrm.est.tau.item( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=tau.item.incr  , theta.k ,
				msteps, mstepconv , tau.item.fixed , prob.rater )
		tau.item <- res$tau.item
		se.tau.item <- res$se.tau.item
		g1  <- abs( tau.item0 - tau.item )
		tau.item.incr <- ifelse( tau.item.incr > g1 , g1 , tau.item.incr )						
		prob.item <- res$prob.item
# cat("est.tau "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;			

		
		# estimate a.item parameter
		if (est.a.item){
		if (iter ==0){	max.b.increment -> a.item.incr }
			res <- .rm.hrm.est.a.item( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=a.item.incr ,theta.k ,
				msteps, mstepconv , prob.rater  , a.item.fixed )		
			a.item <- res$a.item
			se.a.item <- res$se.a.item
			a.item.incr <- abs( a.item0 - a.item )
			prob.item <- res$prob.item
				}

# cat("est.a "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;							
#prob.item <- NULL
				
		# estimate d.rater parameter
		if (est.d.rater!="n"){
			if (iter ==0){	d.rater.incr <- 2 }	
			res <- .rm.hrm.est.d.rater(  c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=d.rater.incr,theta.k ,
					msteps , mstepconv , d.min , d.max , est.d.rater , prob.item  )
			d.rater <- res$d.rater
			se.d.rater <- res$se.d.rater
#			d.rater.incr <- abs( d.rater0 - d.rater )	
#			g1 <- abs( d.rater0 - d.rater )
#			d.rater.incr <- ifelse( d.rater.incr > g1 , g1 , d.rater.incr )									
					}				
# cat("est.d.rater "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;			

					
		# estimate c.rater parameter
		if( est.c.rater!="n" ){	
			if (iter ==0){	max.b.increment -> c.rater.incr }
			res <- .rm.hrm.est.c.rater(  c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=c.rater.incr,theta.k ,
					msteps , mstepconv , est.c.rater , prob.item , c.rater.fixed)
			c.rater <- res$c.rater
			se.c.rater <- res$se.c.rater
			g1 <- abs( c.rater0 - c.rater )
			c.rater.incr <- max(g1)
#			c.rater.incr <- ifelse( c.rater.incr > g1 , g1 , c.rater.incr )						
						}						
 #cat("est.c.rater "); a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 ;			
						
		flush.console()		
		
		
		# update distribution
		res <- rm.smooth.distribution( theta.k , pi.k , est.mean = est.mean , skillspace)
		pi.k <- res$pi.k
		mu <- res$mu
		sigma <- res$sigma
		
		
		dev <- -2*ll
		# convergence criteria
		conv <- max( abs(c.rater-c.rater0) , abs( c.rater-c.rater0) , 
					abs( tau.item0-tau.item) , abs( a.item - a.item0 ) )
		iter <- iter+1
		devchange <- abs( ( dev - dev0 ) / dev0  )
		#****
		# print progress			
		cat( paste( "   Deviance = "  , round( dev , 4 ) , 
			if (iter > 1 ){ " | Deviance change = " } else {""} ,
			if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,"\n",sep="") )
		cat( paste( "    Maximum c.rater parameter change = " , 
				paste( round(max(abs(c.rater0-c.rater)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum d.rater parameter change = " , 
				paste( round(max(abs(d.rater0-d.rater)) ,6) , collapse=" " ) , "\n" , sep=""))								
		cat( paste( "    Maximum tau.item parameter change = " , 
				paste( round(max(abs(tau.item0-tau.item)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum a.item parameter change = " , 
				paste( round(max(abs(a.item0-a.item)) ,6) , collapse=" " ) , "\n" , sep=""))
        cat( paste(" Trait M  = " , round( mu , 3 ) , sep="") , "\n")				
		cat( paste(" Trait SD = " , round( sigma , 3 ) , sep="") , "\n")
		# flush.console()			
				}
				
	# *********
	# arrange OUTPUT
	
	# c parameters
	if ( ! is.null( c.rater.fixed ) ){
		c.rater[ c.rater.fixed[,1:2] ] <- NA
						}
	
	#---
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$VV <- VV
	ic$RR <- RR
	#****
	# skill parameters
	ic$np.skill <- 0
	if ( skillspace == "normal" ){
		ic$np.skill <- 1 + est.mean
						}
	if ( skillspace == "discrete" ){
		ic$np.skill <- length(theta.k) - 1
						}	
	#*****
	# item parameters
	ic$np.item <- VV * max( maxK )
	if ( ! is.null(tau.item.fixed)){
			ic$np.item <- ic$np.item - nrow(tau.item.fixed)
									}
	if ( is.null( a.item.fixed ) ){ 
		ic$np.item <- ic$np.item + est.a.item*(VV-1)
						} else {
		ic$np.item <- ic$np.item + VV - 1
						}
						
	#*****
	# rater parameters
	ic$np.rater <- 0
	if ( est.d.rater=="e" ){ ic$np.rater <- ic$np.rater + 1 }
	if ( est.d.rater=="i" ){ ic$np.rater <- ic$np.rater + VV }	
	if ( est.d.rater=="r" ){ ic$np.rater <- ic$np.rater + RR }	
	if ( est.d.rater=="a" ){ ic$np.rater <- ic$np.rater + I }		
	if ( est.c.rater=="e" ){ 
			ic$np.rater <- ic$np.rater + K 
						}
	if ( est.c.rater=="i" ){ 
			ic$np.rater <- ic$np.rater + sum(maxK)
							}	
	if ( est.c.rater=="r" ){ 
			ic$np.rater <- ic$np.rater + K*RR 
							}		
	if ( est.c.rater=="a" ){ 
			ic$np.rater <- ic$np.rater + K*I 
		if ( ! is.null( c.rater.fixed ) ){
			ic$np.rater <- ic$np.rater - nrow( c.rater.fixed )
						}			
						}		
	ic$np <- ic$np.skill + ic$np.item + ic$np.rater
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	
	#---
	# person
	person <- procdata$person.index
	NP <- nrow(person)
	person$score <- rowSums( dat2 * dat2.resp )
	mkrr <- rep( maxK , RR )
	person$maxscore <- rowSums( dat2.resp * outer( rep(1,NP) , mkrr ) )
	person$EAP <- rowSums( f.qk.yi * outer( rep(1,NP) , theta.k) )
	person$SE.EAP <- sqrt( rowSums( f.qk.yi * outer( rep(1,NP) , theta.k^2) ) - 
			( person$EAP) ^2 )
	EAP.rel <- 1 - mean( person$SE.EAP^2 ) / 
				( mean( person$SE.EAP^2 ) + var( person$EAP ) )

	
	
	#---
	# item
	if (!is.null(tau.item.fixed)){
#		 tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA
#        tau.item[ tau.item.fixed1[,1:2,drop=FALSE] ] <- NA	
		K <- max(maxK)
        I <- nrow(tau.item)
		for (ii in 1:I){
		  if ( maxK[ii] < K ){
			for (kk in seq(maxK[ii]+1,K) ){
				tau.item[ii,kk] <- NA
										}
								}
							}
		se.tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA								
							}

							
    item <- data.frame( "item" = colnames(dat) , 
			"N" = colSums( 1-is.na(dat)) , 
			"M" = colMeans( dat , na.rm=T ) )
	for (kk in 1:K){ item[ , paste0("tau.Cat",kk) ] <- tau.item[,kk] }
    item$a <- a.item
#	item$b <- rowMeans(tau.item)

	# latent mean and standard deviation
	me1 <- rep(NA,VV)
	sd1 <- rep(NA,VV)
	for (ii in 1:VV){
		# ii <- 1
		pii <- prob.item[ii,,]
		qii <- matrix( c(0,Qmatrix[ii,]) , nrow= K+1 , ncol=ncol(pii) )
		me1[ii] <- sum( colSums( qii * pii ) * pi.k )
		sd1[ii] <- sqrt( sum( colSums( qii^2 * pii ) * pi.k ) - me1[ii]^2  )
				}
	item$latM <- me1
	item$latSD <- sd1	
	
	
	obji <- item
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
	cat("*********************************\n")
	cat("Item Parameters\n")
    print( obji )		
	
	#---
	# rater
	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
#	N <- aggregate( N , list( rater.index ) , sum )[,2]
#	M1 <- aggregate( M1 , list( rater.index ) , mean )[,2]	
    rater <- data.frame( "item.rater" = colnames(dat2) , 
			"N" = N , 
			"M" = M1 , 
			"d" = d.rater )
    for (zz in 1:(ncol(c.rater) ) ){
		rater[ , paste0("c_",zz)] <- c.rater[,zz] }
	# transformed c parameters
    for (zz in 1:(ncol(c.rater) ) ){
		# rater[ , paste0("c_",zz,".trans")] <- c.rater[,zz] / d.rater / K 
        rater[ , paste0("c_",zz,".trans")] <- c.rater[,zz] / d.rater	
			}


	rater <- rater[ order( paste( rater$item.rater) ) , ]
	rownames(rater) <- NULL
	rownames(item) <- NULL

	rt1 <- rater
	l1 <- paste(rt1$item.rater)
	l2 <- strsplit( l1 , split="-" )
	rt1$item <- unlist( lapply( l2 , FUN = function(uu){ uu[[1]] } ) )
	rt1$rater <- unlist( lapply( l2 , FUN = function(uu){ uu[[2]] } ) )
	
	#*****
	# dimnames probs
	dimnames(probs)[[1]] <- colnames(dat2)
	#*****			
						
	#*****
	# distribution	
	skill.distribution <- data.frame("theta.k" = theta.k , "pi.k" = pi.k )
	
	#*****
	# labels
	dimnames(prob.item) <- list( colnames(dat) , paste0("Cat" , 0:K) , NULL )
	dimnames(prob.rater) <- list( colnames(dat2) , paste0("Cat" , 0:K) , NULL )	

	
						
	obji <- rater
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
	cat("*********************************\n")
	cat("Rater Parameters\n")		
    print( obji )		
	cat("*********************************\n")
	cat("EAP Reliability = " , round(EAP.rel,3) , "\n")		
	
	s2 <- Sys.time()

	

    res <- list("deviance" = dev , "ic"=ic , "item"=item , "rater"=rater ,
		"person" = person , "EAP.rel"=EAP.rel , 
		"mu"=mu , "sigma"=sigma , 
		"theta.k" = theta.k , "pi.k" = pi.k , "G"=1 , 		
		"tau.item"=tau.item , "se.tau.item"=se.tau.item ,
		"a.item"=a.item , "se.a.item"=se.a.item ,
		"c.rater"=c.rater , "se.c.rater"=se.c.rater , 
		"d.rater"=d.rater , "se.d.rater"=se.d.rater , 
		"f.yi.qk"=f.yi.qk , "f.qk.yi"=f.qk.yi , "probs"=probs ,
		"prob.item" = prob.item ,
		"n.ik"=n.ik , "maxK"=maxK , "pi.k" = pi.k ,
		"procdata" =procdata , "iter"=iter , "theta.k" = theta.k , 
		"Qmatrix" = Qmatrix , 
		"s1"=s1 , "s2"=s2 , "tau.item.fixed"=tau.item.fixed ,
		"rater2" = rt1 , "maxK" = maxK , "skill.distribution"=skill.distribution ,
		"skillspace" = skillspace )
	class(res) <- "rm.sdt"
	return(res)

		}

		
		
################################
# time stamps

# many persons, one item, one rater

# calcprob Time difference of 0 secs
# posterior Time difference of 0 secs
# est.tau Time difference of 0.07799983 secs
# est.a Time difference of 0 secs
# est.c.rater Time difference of 0.04680014 secs
