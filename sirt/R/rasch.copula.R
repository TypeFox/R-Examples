


#------------------------------------------------------------------
# Copula estimation in a raschtype model
##NS export(rasch.copula)
rasch.copula <- function( dat , itemcluster , group = NULL , 
						progress = TRUE , mmliter = 1000 , delta = NULL ,
						theta.k = seq(-4,4,len=21) , alpha1=0 , alpha2=0 , numdiff.parm = .000001 ,
						est.b = seq(1,ncol(dat)) , est.a = rep(1,ncol(dat)) , est.delta = NULL , 
						est.alpha = FALSE , 
						glob.conv = .001 , alpha.conv = .0001 , conv1 = .001 
										){
	###############################################################
	# INPUT:
	# dat	... data frame
	# itemcluster ... vector of integers
	# progress	... display progress?
	# mmliter 	... see rasch.mml
	# delta ... initial delta estimate (is relevant for fixing delta parameters)
	# theta.k... number of grid theta points
	# alpha1, alpha2 ... rasch type parameter
	# numdiff.parm ... numerical differentiation parameter
	# est.b	... which b parameters shall be estimated
	###############################################################
	s1 <- Sys.time()
	# item clusters
    t1 <- table(itemcluster)
	t1b <- as.numeric( names(t1))
	t1c <- t1[ names(t1) != 0 ]
    if ( any( t1c == 1) ){
           stop( "There should be at least two items in an item cluster\n" )
                                    }			
            x1 <- seq( 1 , max(t1b) ) 
            x2 <- sort(setdiff( as.numeric(sort(names(t1))) , 0 ))
            if ( sum( abs( x1-x2) ) > 10^(-10) ){
               stop( "Item cluster identifiers must be recoded to 1, ..., C\n" )
                                 }
	CC <- length(x1) # number of clusters
    # calculation of number of itemclusters	
    if ( progress  ){
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function rasch.copula\n")		
        cat("---------------------------------------------------------------------------------------------------------- \n")
        flush.console()
      }
	I <- ncol(dat)
	if ( is.null( colnames(dat))){
		colnames(dat) <- paste("Item" , 1:I , sep="")
				}
	# groups
	if ( ! is.null( group) ){
		groups <- unique( group )
		G <- length( groups )	} else { G <- 1 }
	GG <- G
	# data preparation (frequencies)
	dat10 <- dat00 <- dat
	dat10[ is.na( dat10) ] <- 9
	patt <- paste("P",dat10[,1] ,sep="")
	for (ii in 2:I){ patt <- paste( patt , dat10[,ii],sep="") }
	pattern <- data.frame( table(patt) )
	colnames(pattern) <- c("pattern" , "freqwgt")
	# patttern in data
#	pattern.in.data <- data.frame(match( patt , pattern$pattern )
	pattern.in.data <- patt
	# calculate frequencies in multiple group case
	if (G > 1 ){ 
		for (gg in 1:G){
		# gg <- 1
		t1 <- table( patt[ group == gg  ] )
		pattern <- merge( pattern , t1 , by.x=1 , by.y = 1 , all=T )
					}
		pattern[ is.na(pattern) ] <- 0	
		colnames(pattern)[-c(1:2)] 	<- paste("freqwgt" , 1:G , sep="")		
				}
	dat0 <- matrix( 0 , nrow=nrow(pattern) , I )
	for (ii in 1:I){ 
			dat0[,ii] <- as.numeric( substring( paste( pattern[,1] ) , ii+1 , ii+1 ) )
				}
	dat0[ dat0 == 9 ] <- NA
	colnames(dat0) <- colnames(dat)
	# define mu and sigma
	mu <- rep(0,G)
	sigma <- rep(1,G)
	if ( G > 1){ 
#			mu[2] <- 1 ; sigma[2] <- 1.2 
				}
	#------------------------
	dat <- dat0
	dat2 <- dat
	dat2.resp <- 1 * ( !is.na( dat2) )			
	dat2[ is.na(dat2) ] <- 0
		# data preparation in case of item clusters
		dat2.ld <- NULL		# dat2.ld is NULL if there are no item clusters
			if ( is.null(delta)){ delta <- stats::runif( CC , .2 , .8 ) }	# initial estimate of delta
			if ( is.null(est.delta)){ est.delta <- 1:CC }
			dp.ld <- as.list( 1:CC )
			# item pattern
			for (cc in 1:CC){	
	#			cc <- 1
				icl.cc <- which( itemcluster == cc )
				dp.ld.cc <- .calc.copula.itemcluster( D = length(icl.cc) )
				dp.ld.cc$items <- icl.cc
				dp.ld.cc$N.items <- NCC <- length(icl.cc)
				dp.ld.cc$itemlabels <- colnames(dat)[icl.cc]
				# item selection for independent items
				m1 <- outer( rep(1,2^NCC) , icl.cc)
				m2 <- dp.ld.cc$patt * m1   + ( 1 - dp.ld.cc$patt ) * ( m1 + I )
				m2 <- matrix( t(m2) ,  nrow=1 , byrow=T )[1,]
				res <- list( "items" = m2 )
				res1 <- rep(NCC , 2^NCC )
				names(res1) <- rownames( dp.ld.cc$patt )
				res$N.Index1 <- res1
				dp.ld.cc$independent <- res
				# item selection for dependent items
				m2 <- ( 1 - dp.ld.cc$patt ) * ( m1 + I )	
				m2 <- matrix( t(m2) ,  nrow=1 , byrow=T )[1,]				
				m2 <- c( m2[ m2> 0 ] , 2*I + 1 )
				res <- list( "items" = m2 )				
				res$N.Index1 <- rowSums( dp.ld.cc$patt  == 0 )
				res$N.Index1[ length(res$N.Index1) ] <- 1
				dp.ld.cc$dependent <- res
				dp.ld[[cc]] <- dp.ld.cc
							}
			# create data frame with item response pattern
			dat2.ld <- matrix(0 , nrow(dat2) , CC )
			dat3.ld <- as.list( 1:CC )
			for (cc in 1:CC){
				# cc <- 1	
				dp.cc <- dp.ld[[cc]]
				dat2.cc <- dat2[ , dp.cc$items ]
				l1 <- "P"
				for ( vv in seq(1,ncol(dat2.cc))){
					l1 <- paste( l1 , dat2.cc[,vv] , sep="")
						}
				dat2.ld[ , cc ] <- match( l1 , rownames( dp.cc$patt ) )
				dat2.ld[ rowSums( dat2.resp[ , dp.cc$items ]	) < length(dp.cc$items) , cc ] <- NA			
				NRR <- nrow( dp.cc$patt )
				dat3.ld.cc <- sapply( seq( 1 , NRR) , FUN = function(aa){ 1*(dat2.ld[,cc] == aa ) } )
				dat3.ld.cc[ is.na(dat3.ld.cc) ] <- 0
				dat3.ld[[cc]] <- dat3.ld.cc
							}					
				# response indicator
				dat2.ld.resp <- 1 - is.na( dat2.ld )
				# set missings in dat2 to some arbitrary category
				dat2.ld[ is.na( dat2.ld ) ] <- 1
				# create dat2 data sets for local independence items
				itemcluster0 <- ind2 <- which( itemcluster == 0 )
				bdat2.li.resp <- dat2.li.resp <- bdat2.li <- dat2.li <- NULL
				
				if ( length(ind2) > 0 ){ 
						dat2.li <- dat2[ , ind2 ]
						dat2.li.resp <- dat2.resp[ , ind2 ]
						bdat2.li <- cbind( dat2.li , 1 - dat2.li )
						bdat2.li.resp <- cbind( dat2.li.resp , dat2.li.resp )
										}
		#########################################################################
	#--------------------------------------------------
	# initial estimate of item difficulty
#	b <- rasch.pairwise( dat , progress = FALSE)$item$itemdiff
	b <- - stats::qlogis( ( colMeans( dat00 , na.rm=T ) + .005 ) / 1.01 )
	I <- ncol(dat2)
	# initial estimate of (mean) item discrimination
	a <- rep( 1 , I )
	# density weights
	wgt.theta <- stats::dnorm(theta.k)
	wgt.theta <- wgt.theta / sum( wgt.theta )
    if ( G > 1){
		wgt.theta <- matrix(0 , length(theta.k) , G )
		for ( gg in 1:G){
			wgt.theta[,gg] <- stats::dnorm( theta.k , mean = mu[gg] , sd = sigma[gg] )
			wgt.theta[,gg] <- wgt.theta[,gg] / sum( wgt.theta[,gg] )			
						}
				}
	iter <- 0
	#**********************************
	# BEGIN MARGINAL MAXIMUM LIKELIHOOD ESTIMATION
	dev <- 1 ; par.change <- dev.change <- 1000 
	maxalphachange <- 0
	while ( ( dev.change > glob.conv | par.change > conv1 | maxalphachange > alpha.conv ) & iter < mmliter ){
		
		cat( paste(rep("-" , 70), collapse="") , "\n")
		k1 <- floor( log10(iter+1) )
		x1 <- "        |" 
		x1 <- substring( x1 , k1+1 )
		s1c <- Sys.time()
		cat( paste( paste( "MML EM Iter." , iter + 1 ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
						s1c  , "  " ,
						if ( iter > 0 ){ paste( round(difftime(s1c ,s1b , units='secs' ),4) , "secs" ) } , 
						"\n" ,sep="") ) # 
		s1b <- Sys.time()
		h <- numdiff.parm 
		dev0 <- dev
		#************************************
		# estimation of b parameters
		b0 <- b
		# identify different b parameter groups
		bG <- setdiff( unique( est.b ) , 0 )
		prbar <- seq( 1 , 10 , len = length(bG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of b:     |")
		for (bb in bG){
			est.bb <- 1 * (est.b == bb )
			if (bb == 1 ){ 
				rescop <- .ll.rasch.copula2( theta.k , b0 , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG  )
				ll0 <- rescop$ll
							}							
			ll1 <- .update.ll.rasch.copula( theta.k , b0 + h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG )$ll
			ll2 <- .update.ll.rasch.copula( theta.k , b0 - h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG )$ll

			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
			b.change <- - d1 / d2
			b.change <- ifelse( abs( b.change ) > .5 , .5*sign(b.change) , b.change )              
			b.change <- b.change * est.bb
			b <- b + b.change
#			cat( bb , " ") ; 
			cat( paste( rep( "-" , prbar[bb]), collapse="") )
			flush.console()				
					}
        a1b <- max( abs( b - b0 ) )
		cat("|     max. parm. change" , round( a1b , 5),"\n")
		#******************************************************************************
		# estimation of a parameters
		a0 <- a
		# identify different a parameter groups
		aG <- setdiff( unique( est.a ) , 0 )
		prbar <- seq( 1 , 10 , len = length(aG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of a:     |")
		for (aa in aG){
			est.aa <- 1 * (est.a == aa )
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll0 <- rescop$ll
			ll1 <- .update.ll.rasch.copula( theta.k , b, alpha1 , alpha2 , a + h*est.aa , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG )$ll
			ll2 <- .update.ll.rasch.copula( theta.k , b , alpha1 , alpha2 , a - h*est.aa , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			a.change <- - d1 / d2
			a.change <- ifelse( abs( a.change ) > .3 , .3*sign(a.change) , a.change )              
			a.change <- a.change * est.aa
			a <- a + a.change
			a[ a < 0 ] <- .01			
#			cat( aa , " ") ; 
			cat( paste( rep( "-" , prbar[aa]), collapse="") )
			flush.console()
							}
		if ( length(aG) < 2 ){ cat( paste( rep( "-" , 10 - length(aG) ), collapse="") ) }
		a1a <- max( abs( a - a0 ) )
		cat("|     max. parm. change" , round( a1a , 5),"\n")
		#******************************************************************************
		# estimation of delta parameters
		delta0 <- delta
		dG <- setdiff( unique( est.delta ) , 0 )
		prbar <- seq( 1 , 10 , len = length(dG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		# identify different a parameter groups
		cat(" Estimation of delta: |")
		for (cc in dG){
			est.cc <- 1 * ( est.delta == cc )
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , 
							delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern, GG)
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula( theta.k , b, alpha1 , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta + h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG)$ll
			ll2 <- .update.ll.rasch.copula( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta - h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			delta.change <- - d1 / d2
			delta.change <- ifelse( abs( delta.change ) > .2 , .2*sign(delta.change) , delta.change )              
			delta.change <- delta.change * est.cc
			delta <- delta + delta.change
			delta[ delta <= 0 ] <- 2*numdiff.parm			
			delta[ delta >= 1 ] <- 1 - 2*numdiff.parm	
# adapt delta boundary estimate
#if ( adapt.delta.boundary & ( iter %% 5 == 0 ) ){
# print("TRUE")
#	delta <- ifelse ( abs( delta - 2*numdiff.parm ) < 10^(-15) , .20 , delta )
#				}
			
			cat( paste( rep( "-" , prbar[cc]), collapse="") )
			flush.console()
							}
		if ( length(dG) < 2 ){ cat( paste( rep( "-" , 10 - length(dG) ), collapse="") ) }
        a1d <- max( abs( delta - delta0 ) )
		cat("|     max. parm. change" , round( a1d , 5),"\n")
		#******************************************************************************
		# estimation of alpha parameters
		alpha10 <- alpha1
		alpha20 <- alpha2
		prbar <- 5
		cat(" Estimation of alpha: |")		
		# alpha1
		if (est.alpha){
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG)
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula( theta.k , b, alpha1 + h , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG)$ll
			ll2 <- .update.ll.rasch.copula( theta.k , b , alpha1 - h , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k1 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha1 <- alpha1 + alpha.change
			}
			cat( paste( rep( "-" , prbar), collapse="") )
			flush.console()		
		# alpha2
		if (est.alpha){
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG )
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula( theta.k , b, alpha1 , alpha2+h , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG)$ll
			ll2 <- .update.ll.rasch.copula( theta.k , b , alpha1 , alpha2 -h, a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k2 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha2 <- alpha2 + alpha.change
			}
			cat( paste( rep( "-" , prbar), collapse="") )
			flush.console()		
        a1k <- max( abs( c( alpha1 - alpha10, alpha2 - alpha20 )) )
		cat("|     max. parm. change" , round( a1k , 5),"\n")
		#******************************************************************************
		# estimation of mu parameters
		a1m <- 0
		if (G>1){
		mu0 <- mu
		# identify different a parameter groups
#		h <- 10 * numdiff.parm
#		dG <- setdiff( unique( est.delta ) , 0 )
#		prbar <- seq( 1 , 10 , len = length(dG) )
#		prbar <- floor( prbar )
#		prbar <- c( prbar[1] , diff(prbar) )
#		muG <- GG-1
		muG <- 1:(GG-1)
		prbar <- seq( 1 , 10 , len = length(muG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of mu:    |")
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll0 <- rescop$ll
			# mu + h		
			w1 <- wgt.theta
			w2 <- stats::dnorm( theta.k , mean = mu[gg] + h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll1 <- rescop$ll
			# mu - h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg] - h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll2 <- rescop$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			mu.change <- - d1 / d2
			mu.change <- ifelse( abs( mu.change ) > .3 , .3*sign(mu.change) , mu.change )      		
			mu.change <- mu.change * ( ( 1:G ) == gg )
			mu <- mu + mu.change
			w2 <- stats::dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] )
			wgt.theta[,gg] <- w2 / sum(w2)
#			cat( aa , " ") ; 
			cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
			flush.console()
							}
		if ( length(muG) < 2 ){ cat( paste( rep( "-" , 10 - length(muG) ), collapse="") ) }
		a1m <- max( abs( mu - mu0 ) )
		cat("|     max. parm. change" , round( a1m , 5),"\n")
				}			# end mu
		######################################################################
		#******************************************************************************
		# estimation of sigma parameters
		a1s <- 0
		if (G>1){
		sigma0 <- sigma
		# identify different a parameter groups
#		h <- 10 * numdiff.parm
		sigmaG <- seq(1,GG-1)
		prbar <- seq( 1 , 10 , len = length(sigmaG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of sigma: |")
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll0 <- rescop$ll
			# sigma + h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] +h)
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll1 <- rescop$ll
			# sigma - h		
			w1 <- wgt.theta
			w2 <- stats::dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg]-h )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG )
			ll2 <- rescop$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			sigma.change <- - d1 / d2
			sigma.change <- ifelse( abs( sigma.change ) > .3 , .3*sign(sigma.change) , sigma.change )      		
			sigma.change <- sigma.change * ( ( 1:G ) == gg )
			sigma <- sigma + sigma.change
			w2 <- stats::dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] )
			wgt.theta[,gg] <- w2 / sum(w2)
#			cat( aa , " ") ; 
			cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
			flush.console()
							}
		if ( length(sigmaG) < 2 ){ cat( paste( rep( "-" , 10 - length(sigmaG) ), collapse="") ) }
		a1s <- max( abs( sigma - sigma0 ) )
		cat("|     max. parm. change" , round( a1s , 5),"\n")
				}			# end sigma
		######################################################################		
		
		
		
		#**********************************************************************************
		# convergence display 
		a1 <- stats::aggregate( b , list( est.b) , mean )
		cat("   b parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )
		a1 <- stats::aggregate( a , list( est.a) , mean )
		cat("   a parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )
		cat("   delta parameters: " , paste( round( delta , 3 ) , collapse= " " ) , "\n" )
		cat("   alpha parameters: " , paste( round( c(alpha1 , alpha2) , 3 ) , collapse= " " ) , "\n" )		
		cat("   mu parameters: " , paste( round( mu , 3 ) , collapse= " " ) , "\n" )		
		cat("   sigma parameters: " , paste( round( sigma , 3 ) , collapse= " " ) , "\n" )
		#******************************************************************************
		iter <- iter + 1 
		dev <- -2*ll0
        dev.change <- abs( ( dev - dev0)/ dev0 )
        par.change <- max( a1a , a1b , a1d , a1k , a1m , a1s)
		cat( "Deviance = "  ,   round( dev , 5 ) , "| max. parm. change = " , 
                                        round( par.change , 6 ) ,  " \n"   )  
		if ( ( dev > dev0 ) & ( iter > 4 ) ){ cat("   Deviance has increased! Convergence Problems?\n") }
		flush.console()
#stop("herer")
			}
	# end MML iterations
	#**********************************************************************************
	# Standard error estimation (This is a TO DO!)
	iterend <- iter
#	iter <- 1
#	cat( paste(rep("-" , 70), collapse="") , "\n")
#	k1 <- floor( log10(iter+1) )
#	x1 <- " |" 
#	x1 <- substring( x1 , k1+1 )	
#	cat( paste( paste( "Standard errors (SE's)"  ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
#					Sys.time() , "\n" ,sep="") ) # 
#	h <- numdiff.parm 
#	dfr <- NULL
		#************************************
		# standard errors b
		# identify different b parameter groups
#		bG <- setdiff( unique( est.b ) , 0 )
#		prbar <- seq( 1 , 10 , len = length(bG) )
#		prbar <- floor( prbar )
#		prbar <- c( prbar[1] , diff(prbar) )
#		cat(" SE's of b:            |")	
#		for (bb in bG){
#			est.bb <- est.b * (est.b == bb )
#			if (bb == 1 ){ 
#				rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp )
#				ll0 <- rescop$ll
#							}							
#
#			ll1 <- .update.ll.rasch.copula( theta.k , b + h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp , rescop , itemcluster )
							
							
#			ll2 <- .update.ll.rasch.copula( theta.k , b - h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp , rescop , itemcluster )
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
#			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
#			d2 <- sqrt( 1 / sum( - d2 ) )
#			dfr.pp <- data.frame( "parameter" = "b" , "pargropup" = bb , 
#							"est" = b[ bb ] , "se" = d2 )
#			dfr <- rbind( dfr , dfr.pp )
#			cat( paste( rep( "-" , prbar[bb]), collapse="") )
#			flush.console()
#							}
#		cat("|\n")
#    print(dfr)

	#********************************************************
	# information criteria
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat00) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
        ic[[ "np" ]] <- length(bG) + length(aG) + length(dG) + 2*est.alpha
        # AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np	
	#**********************************************************************************
	# results item parameters			
	item <- data.frame( "item" = colnames(dat) , 
				"N" = colSums(!is.na(dat00)) , 
				"p" = colMeans( dat00 , na.rm=T ), 
				"b" = b , "est.b"= est.b , "a" = a , "est.a" = est.a )
	item$thresh <- item$a * item$b
	# add results dependency parameter for item clusters
	item$itemcluster <- itemcluster
	item$delta <- 0
    cat("---------------------------------------------------------------------------------------------------------- \n")
	for (cc in 1:CC){
		# cc <- 1
		dcc <- dp.ld[[cc]]
		item[ dcc$items , "delta"] <- delta[cc]
				}
	item$est.delta <- 0
	i1 <- which( item$itemcluster > 0 )
	icl <- item$itemcluster[ i1 ]		
	item$est.delta[ i1 ] <- est.delta[ icl ]	
	# print item summary
	cat("Parameter summary\n")
	.pr( item , digits=3 )		# print item statistics
#	print(round(item,3))
        # computation time
        s2 <- Sys.time()
        if (progress){ 
                cat("---------------------------------------------------------------------------------------------------------- \n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("---------------------------------------------------------------------------------------------------------- \n")
                    }      
    # collect results

#	 dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp 
	datalist <- list( pattern.in.data = pattern.in.data , dat0 = dat0 ,
					dat2 = dat2 , dat2.resp = dat2.resp , dat2.li = dat2.li , 
					dat2.ld = dat2.ld , dat2.li.resp=dat2.li.resp , 
					dat2.ld.resp = dat2.ld.resp , dp.ld = dp.ld , CC = CC ,
					bdat2.li = bdat2.li , bdat2.li.resp = bdat2.li.resp , 
					itemcluster0 = itemcluster0 , dat3.ld = dat3.ld		
							)						
	res <- list( "N.itemclusters" = CC , "item" = item , "iter" = iterend , "dev" = dev ,
					"delta" = delta , "b" = b , "a" = a , "mu" = mu , "sigma" = sigma , 
					"alpha1"=alpha1 , "alpha2"=alpha2 , "ic" = ic , "theta.k" = theta.k , "deviance" = dev ,
					"pattern" = pattern	 , "datalist" = datalist				)	
	class(res) <- "rasch.copula"
	return(res)
		}
#----------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
		ind.b <- which( b != rescopula$b )
		ind.a <- which( a != rescopula$a )
		ind.delta <- which( delta != rescopula$delta )
		ind.alpha1 <- ( alpha1 != rescopula$alpha1 )	+ ( alpha2 != rescopula$alpha2 )
		if (ind.alpha1 > 0){ ind.alpha <- seq(1 , ncol(dat2.ld) ) } else { ind.alpha <- NULL }
		itemset <- union( ind.b , ind.a )
		itemset <- union( itemset , ind.alpha )	
		# update term local independence
		li.update <- 1 * ( sum( itemcluster0 %in% itemset ) > 0 )
		# update terms item dependence parameters
		ld.update <- sapply( 1:CC , FUN = function(cc){ 
				g1 <- intersect( which( itemcluster == cc )  , itemset )
				if ( length(g1)){ v1 <- cc } else { v1 <- NULL }
				v1
					} )
		ld.update <- unique( union( ind.delta , unlist( ld.update) ) )
		#------------------------------------------	
		ndat2 <- nrow(dat2.ld)
		ntheta <- length(theta.k)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k.tt <- pjk.theta.k[ rep( seq(1,ntheta) , each = ndat2 ) , ]
		# define indices for persons 
		ind.dat2 <- rep( seq(1,ndat2) , ntheta )
		g1 <- rescopula$g1
		M1 <- rep(1,ndat2)
		pqjk.theta.k.tt <- cbind( pjk.theta.k.tt , 1 - pjk.theta.k.tt , 1 )
			#*************************************************************
			if ( ( length(itemcluster0) > 0 ) & ( li.update > 0 ) ){
				pqjk.theta.k.tt0 <- pqjk.theta.k.tt[ , c( itemcluster0 , itemcluster0+I) ]
				# likelihood for independent items at theta tt
				ll.tt <- ( pqjk.theta.k.tt0 ^ bdat2.li[ind.dat2,] )^bdat2.li.resp[ind.dat2,]
#				g1 <- rowProds( ll.tt ) 		} else { g1 <- 1 }
				g1[,1] <- rowProds2( ll.tt ) 		} # else { g1[,1] <- 1 }
				#*************************************************************
				# likelihood for dependent items
			for (cc in ld.update){
				# cc <- 1
				dat3.ld.cc <- dat3.ld[[cc]] 
				dat3.ld.cc <- dat3.ld.cc[ ind.dat2 , ]
				dp.ld.cc <- dp.ld[[cc]]
				m1.cc <- pqjk.theta.k.tt[ , dp.ld.cc$independent$items ]
				v1.cc <- dp.ld.cc$independent$N.Index1
				# product under independence				
#				Fpjk.cc <- .rowProds.bundle( m1 = m1.cc , v1 = v1.cc)
				Fpjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				#print( dim(m1.cc))
				# evaluate likelihood
				g1.tt <- ( rowSums(Fpjk.cc * dat3.ld.cc ) )^dat2.ld.resp[ind.dat2 ,cc]
				# product under dependence	
				m1.cc <-  pqjk.theta.k.tt[ , dp.ld.cc$dependent$items ] 
				v1.cc <- dp.ld.cc$dependent$N.Index1
#				F0pjk.cc <- .rowMins.bundle( m1 = m1.cc , v1 = v1.cc)
				F0pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
				F1pjk.cc <- F0pjk.cc %*% t(dp.ld.cc$calc)
				g2.tt <- ( rowSums(F1pjk.cc * dat3.ld.cc) )^dat2.ld.resp[ind.dat2,cc]
				g3.tt <- ( 1 - delta[cc] ) * g1.tt + delta[cc] * g2.tt
				g1[,cc+1] <- g3.tt
#				g1 <- g1*g3.tt		# end product for item cluster cc
							}
				res <- g1		
				g1 <- rowProds2( g1 )
				ll.theta.k <- matrix( g1 , ncol= ntheta )
				# casewise: integrate over the product
				# print(ll.theta.k)
				# print( outer( M1 , wgt.theta ) )
#				lli <- log( rowSums(  ll.theta.k * outer( M1 , wgt.theta )  ) ) 
#				ll1 <- sum( pattern$freqwgt * lli )
				# one group
				if (G == 1){ 
					# casewise: integrate over the product evaluate at grid of trait density
					lli <- log( rowSums(  ll.theta.k * outer( M1 , wgt.theta )  ) )
					ll1 <- sum( lli * pattern$freqwgt)
							}
				# mutiple groups
				if (G>1){
					ll1 <- 0
					lli <- as.list( rep(1,G) )
					for (gg in 1:G){
						lli[[gg]] <- log( rowSums(  ll.theta.k * outer( M1 , wgt.theta[,gg] )  ) )
						ll1 <- ll1 + sum( lli[[gg]] * pattern[,gg+2] )
									}
							}				
			
#				ll1 <- sum( log( rowSums(  ll.theta.k * outer( M1 , wgt.theta )  ) ) )
				res <- list( "ll" = ll1 , "lli" = lli )
				return(res)
				}
#----------------------------------------------------------------------------------------------				




#----------------------------------------------------------------------------------------------
.ll.rasch.copula2 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , pattern , GG , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		G <- GG
		ndat2 <- nrow(dat2.ld)
		ntheta <- length(theta.k)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k.tt <- pjk.theta.k[ rep( seq(1,ntheta) , each = ndat2 ) , ]
		# define indices for persons 
		ind.dat2 <- rep( seq(1,ndat2) , ntheta )
		g1 <- matrix( 0 , nrow(pjk.theta.k.tt) , CC + 1 )
		M1 <- rep(1,ndat2)
		pqjk.theta.k.tt <- cbind( pjk.theta.k.tt , 1 - pjk.theta.k.tt , 1 )
			if ( length(itemcluster0) > 0 ){
				pqjk.theta.k.tt0 <- pqjk.theta.k.tt[ , c( itemcluster0 , itemcluster0+I) ]
				# likelihood for independent items at theta tt
				ll.tt <- ( pqjk.theta.k.tt0 ^ bdat2.li[ind.dat2,] )^bdat2.li.resp[ind.dat2,]
#				g1 <- rowProds( ll.tt ) 		} else { g1 <- 1 }
				g1[,1] <- rowProds2( ll.tt ) 		} else { g1[,1] <- 1 }
				# likelihood for dependent items
			for (cc in 1:CC){
				# cc <- 1
				dat3.ld.cc <- dat3.ld[[cc]] 
				dat3.ld.cc <- dat3.ld.cc[ ind.dat2 , ]
				dp.ld.cc <- dp.ld[[cc]]
				m1.cc <- pqjk.theta.k.tt[ , dp.ld.cc$independent$items ]
				v1.cc <- dp.ld.cc$independent$N.Index1
				# product under independence				
#				Fpjk.cc <- .rowProds.bundle( m1 = m1.cc , v1 = v1.cc)
				Fpjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# evaluate likelihood
				g1.tt <- ( rowSums(Fpjk.cc * dat3.ld.cc ) )^dat2.ld.resp[ind.dat2 ,cc]
				# product under dependence	
				m1.cc <-  pqjk.theta.k.tt[ , dp.ld.cc$dependent$items ] 
				v1.cc <- dp.ld.cc$dependent$N.Index1
#				F0pjk.cc <- .rowMins.bundle( m1 = m1.cc , v1 = v1.cc)
				F0pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
				F1pjk.cc <- F0pjk.cc %*% t(dp.ld.cc$calc)
				g2.tt <- ( rowSums(F1pjk.cc * dat3.ld.cc) )^dat2.ld.resp[ind.dat2,cc]
				g3.tt <- ( 1 - delta[cc] ) * g1.tt + delta[cc] * g2.tt
				g1[,cc+1] <- g3.tt
#				g1 <- g1*g3.tt		# end product for item cluster cc
							}
				res <- g1	
				g1 <- rowProds2( g1 )
				ll.theta.k <- matrix( g1 , ncol= ntheta )
				# one group
				if (G == 1){ 
					# casewise: integrate over the product evaluate at grid of trait density
					lli <- log( rowSums(  ll.theta.k * outer( M1 , wgt.theta )  ) )
					ll1 <- sum( lli * pattern$freqwgt)
							}
				# mutiple groups
				if (G>1){
					ll1 <- 0
					lli <- as.list( rep(1,G) )
					for (gg in 1:G){
						lli[[gg]] <- log( rowSums(  ll.theta.k * outer( M1 , wgt.theta[,gg] )  ) )
						ll1 <- ll1 + sum( lli[[gg]] * pattern[,gg+2] )
									}
							}
				res <- list( "ll"=ll1 , "g1"=res , "b" = b , "a" = a , "delta"=delta ,
									"alpha1" = alpha1 , "alpha2" = alpha2 , "lli" = lli )
				}
#----------------------------------------------------------------------------------------------				



