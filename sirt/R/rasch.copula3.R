


#------------------------------------------------------------------
# Copula estimation in a rasch type model
rasch.copula3 <- function( dat , itemcluster , dims=NULL ,
						copula.type = "bound.mixt" ,
						progress = TRUE , mmliter = 1000 , delta = NULL ,
						theta.k = seq(-4,4,len=21) , alpha1=0 , alpha2=0 , 
						numdiff.parm = .000001 ,
						est.b = seq(1,ncol(dat)) , est.a = rep(1,ncol(dat)) , est.delta = NULL , 
						b.init = NULL , a.init = NULL , 
						est.alpha = FALSE , 
						glob.conv = .0001 , alpha.conv = .0001 , conv1 = .001 ,
						dev.crit = .2 , rho.init = .5 , increment.factor = 1.01
#						pattern.off = FALSE
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
	group <- NULL	
	N <- nrow(dat)
	if (is.null(dims) ){
		dims <- rep( 1 , ncol(dat) )
					}	
	
	# arrange item clusters item clusters
    t1 <- table(itemcluster)
	itemcluster[ which( itemcluster %in% names(t1)[ t1 == 1	 ] ) ] <- 0
	itemcluster <- match( itemcluster , unique( sort(itemcluster[itemcluster!=0]) ) )
	itemcluster[ is.na( itemcluster ) ] <- 0
	t1 <- table(itemcluster)	
	t1c <- t1[ names(t1) != 0 ]
	t1b <- as.numeric( names(t1))
	t1c <- t1[ names(t1) != 0 ]
	theta.k <- as.matrix(theta.k)
	I <- ncol(dat)
	
	
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
        cat("-----------------------------------------------------------------\n")
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function 'rasch.copula3'\n")
        cat("-----------------------------------------------------------------\n")
        flush.console()
      }
	# dimensions
    dims <- match( dims , unique( dims ) )
	D <- length(unique(dims))
	theta.k <- as.matrix( expand.grid( as.data.frame( matrix( rep( theta.k, D) , ncol = D ) ) ) )
#	D <- ncol(theta.k)	
	
    G <- 1	# function can only be applied for single groups	
	  
	  
	 # arrange copula types
	if ( length( copula.type ) == 1 ){ copula.type <- rep( copula.type , CC ) }
	  
	I <- ncol(dat)
	if ( is.null( colnames(dat))){
		colnames(dat) <- paste("Item" , 1:I , sep="")
				}
	# remove cases where all responses are missings
	ind <- which( rowMeans( is.na(dat) ) < 1 )
	N0 <- nrow(dat)
	dat <- dat[ ind , ]
	Nmiss <- N0 - nrow(dat)
	if (Nmiss > 0){ 
		cat(paste("Removed ", Nmiss , " cases with only missing reponses\n",sep=""))
				}
	# groups
	if ( ! is.null( group) ){
		cat("Multiple groups are ignored!\n")
		cat("This option is only implemented in 'rasch.copula'\n")
#		groups <- unique( group )		
#		G <- length( groups )	
		G <- 1 ; group <- NULL
			} else { G <- 1 }
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
	# in multidimensional case
	if (D > 0){	
		mu <- rep(0,D) 
		rho <- rho.init
		sigma <- rho + (1-rho) *diag(1,D)
			}
	
				
	#------------------------
	dat <- dat0
	dat2 <- dat
	dat2.resp <- 1 * ( !is.na( dat2) )			
	dat2[ is.na(dat2) ] <- 0
		# data preparation in case of item clusters
		dat2.ld <- NULL		# dat2.ld is NULL if there are no item clusters
			if ( is.null(delta)){ 
					delta <- stats::runif( CC , .3 , .7 ) 
					delta <- ifelse( copula.type == "frank" , 1.3 , delta )
#					delta <- ifelse( copula.type == "cook.johnson" , 1 , delta )										
							}	# initial estimate of delta
			if ( is.null(est.delta)){ est.delta <- 1:CC }
			dp.ld <- as.list( 1:CC )
			# item pattern
			for (cc in 1:CC){	
	#			cc <- 1
				icl.cc <- which( itemcluster == cc )
				dp.ld.cc <- .calc.copula.itemcluster2( D = length(icl.cc) )
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
						dat2.li <- dat2[ , ind2 , drop=FALSE]
						dat2.li.resp <- dat2.resp[ , ind2 , drop=FALSE]
						bdat2.li <- cbind( dat2.li , 1 - dat2.li )
						bdat2.li.resp <- cbind( dat2.li.resp , dat2.li.resp )
										}
		# descriptives itemcluster
		Ncat.ld <- max( unlist( lapply( dp.ld , FUN = function(ll){ nrow(ll$patt) } ) ))
				
		# design table for estimating item difficulties
		b.design <- NULL
		if (length(itemcluster0) > 0){
			b.design <- data.frame( "itemcluster"=0 , 
					"item" = itemcluster0 , "b.indexgroup" = 1 )
									}
		for (cc in 1:CC){
		#	cc <- 1
			g.cc <- ( dp.ld[[cc]] )$items
			bb <- data.frame("itemcluster"= cc , 
						"item" = g.cc , "b.indexgroup" = seq( 1 , length(g.cc) ) )
			b.design <- rbind( b.design , bb )
					}
		b.design <- b.design[ order( b.design$item ) , ]
		b.design$est.b <- est.b
		t1 <- table( setdiff( est.b , 0 ) )
		if ( max(t1) > 1 ){	b.design$b.indexgroup <- 1:I }
		if ( sum( est.b == 0 ) > 0 ){ 
				b.design[ est.b == 0 , "b.indexgroup"  ] <- b.design
							}
		b1 <- setdiff( 1:max(b.design$b.indexgroup) , 
						setdiff( b.design$b.indexgroup , 0 ) )
		if ( length(b1) > 0 ){ 	b.design$b.indexgroup <- 1:I }						
		
	#########################################################################
	#--------------------------------------------------
	# initial estimate of item difficulty
#	b <- rasch.pairwise( dat , progress = FALSE)$item$itemdiff
	I <- ncol(dat2)
	if ( is.null( b.init) ){
		b <- - qlogis( ( colMeans( dat00 , na.rm=T ) + .005 ) / 1.01 )
				} else { b <- b.init }
	# initial estimate of (mean) item discrimination
	if ( is.null(a.init) ){ a <- rep( 1 , I ) } else { a <- a.init }
	# density weights
#	wgt.theta <- dnorm(theta.k)
	wgt.theta <- mvtnorm::dmvnorm( theta.k , mean=mu , sigma=sigma)
	wgt.theta <- wgt.theta / sum( wgt.theta )	
	
    if ( G > 1){
		wgt.theta <- matrix(0 , length(theta.k) , G )
		for ( gg in 1:G){
			wgt.theta[,gg] <- stats::dnorm( theta.k , mean = mu[gg] , sd = sigma[gg] )
			wgt.theta[,gg] <- wgt.theta[,gg] / sum( wgt.theta[,gg] )			
						}
				}
	iter <- 0
	# preliminaries
	M2 <- outer( rep(1,nrow(pattern)), wgt.theta )	
	# define Qmatrix
	Qmatrix <- matrix( 0 , I , D )

	Qmatrix[ cbind(1:I , dims ) ] <- 1
	# define mu.fixed
	mu.fixed <- cbind( 1:D , rep(0,D) )
	# maximum increments
	hstep_b <- .6
	hstep_delta <- .2


	
	#************************************************************************
	#************************************************************************
	# BEGIN MARGINAL MAXIMUM LIKELIHOOD ESTIMATION
	dev <- 1 ; absdev.change <- par.change <- dev.change <- 1000 
	res.posterior <- NULL
	maxalphachange <- 0
	while ( ( ( absdev.change > dev.crit | dev.change > glob.conv | 
					par.change > conv1 | maxalphachange > alpha.conv ) & iter < mmliter ) ){

				
# zz0 <- Sys.time()		
	if (progress){ 
		cat( paste(rep("-" , 70), collapse="") , "\n")
		k1 <- floor( log10(iter+1) )
		x1 <- "        |" 
		x1 <- substring( x1 , k1+1 )
		s1c <- Sys.time()
		cat( paste( paste( "MML EM Iter." , iter + 1 ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
						s1c  , "  " ,
						if ( iter > 0 ){ paste( round(difftime(s1c ,s1b , units='secs' ),4) , "secs" ) } , 
						"\n" ,sep="") ) # 
						}
		s1b <- Sys.time()
		h <- numdiff.parm 
		dev0 <- dev
		#************************************
		# estimation of b parameters
		b0 <- b
		# identify different b parameter groups
#		bG <- setdiff( unique( est.b ) , 0 )
		bG <- setdiff( unique( b.design$b.indexgroup ) , 0 )		
		prbar <- seq( 1 , 10 , len = length(bG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		if (progress){ cat(" Estimation of b:     |")		 }
		for (bb in bG){
#vv0 <- Sys.time()		
			est.bb <- 1 * ( b.design$b.indexgroup == bb )
			b.design.bb <- b.design[ b.design$b.indexgroup == bb , ]
			if (bb == 1 ){ 
				rescop <- .ll.rasch.copula320( theta.k , b0 , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta ,
							wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG , copula.type , 
							Ncat.ld	, Qmatrix=Qmatrix )
				res.posterior <- rescop	
							}
			# is this really necessary?
#			wgt.theta <- rescop$pik
						
			rest1 <- .update.ll.rasch.copula321( theta.k , b0 + h*est.bb , alpha1 , alpha2 , 
				a , dat2.li , itemcluster0 , 
					CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
					bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type,
					Qmatrix=Qmatrix)
## continue here the code rewriting!!!

					
			rest2 <- .update.ll.rasch.copula321( theta.k , b0 - h*est.bb , alpha1 , alpha2 , 
					a , dat2.li , itemcluster0 , 
					CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
					bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type,
					Qmatrix=Qmatrix)
			ll0 <- ll1 <- ll2 <- rep(0,I)
			# numerical derivatives independent items
			if ( rescop$calc.ind ){
			for (jj in 1:(length(itemcluster0) ) ){
				ll0[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rescop$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rescop$pjk.theta.k0[,jj] ) )
				ll1[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rest1$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rest1$pjk.theta.k0[,jj] ) )
				ll2[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rest2$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rest2$pjk.theta.k0[,jj] ) )								
										} }
			itemcluster1 <- b.design.bb[ b.design.bb$itemcluster > 0 , "item"]
			rjkCC <- rescop$rjkCC
			for (jj in 1:(length(itemcluster1) ) ){
					cc <- b.design.bb[ b.design.bb$item == itemcluster1[jj] , "itemcluster" ]
					ll0[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rescop$pjk.theta.kCC[[cc]] ) )
					ll1[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rest1$pjk.theta.kCC[[cc]] ) )
					ll2[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rest2$pjk.theta.kCC[[cc]] ) )			
								}
			a1 <-stats::aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.b) , sum , na.rm=TRUE)					
			ll0 <- a1[,2]
			ll1 <- a1[,3]
			ll2 <- a1[,4]			
			b.change <- nr.numdiff( ll0=ll0 , ll1=ll1 , ll2=ll2 , h=h )	
#			hstep <- .5^( log(iter) )
			if (bb == bG[1] ){
				hstep_b <- hstep <- hstep_b * ( 1 / increment.factor )
							}
			b.change <- ifelse( abs( b.change ) > hstep , hstep*sign(b.change) , b.change )            			
			b.change <- b.change[ match( est.b , a1[,1] ) ]		
			b <- b + b.change
			if (progress){ 
				cat( paste( rep( "-" , prbar[bb]), collapse="") ) 
				flush.console()	
						}
# cat("end b item") ; vv1 <- Sys.time(); print(vv1-zz0) ; vv0 <- vv1				
					}
        a1b <- max( abs( b - b0 ) )
		if (progress){
			cat("|     max. parm. change" , round( a1b , 5),"\n")
					}
		
		wm1 <- sum( theta.k * rescop$pik )
		wsd <- sqrt( sum( ( theta.k - wm1 )^2 * rescop$pik ) )	
# cat("end b") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
		#******************************************************************************
		# estimation of a parameters
		a0 <- a
		# identify different a parameter groups
		aG <- setdiff( unique( est.a ) , 0 )
		prbar <- seq( 1 , 10 , len = length(aG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		if (progress){ cat(" Estimation of a:     |") }
		for (aa in aG){
			est.aa <- 1 * (est.a == aa )	
			rescop <- .ll.rasch.copula320( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , 
							delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type , Ncat.ld,
							Qmatrix=Qmatrix)
			ll0 <- rescop$ll
			ll1 <- .update.ll.rasch.copula320( theta.k , b, alpha1 , alpha2 , a + h*est.aa , 
						dat2.li , itemcluster0 , 
						CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
						bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type ,
						Qmatrix=Qmatrix)$ll
			ll2 <- .update.ll.rasch.copula320( theta.k , b , alpha1 , alpha2 , a - h*est.aa , 
						dat2.li , itemcluster0 , CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I ,	bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type, Qmatrix=Qmatrix)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )  
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2	
 			if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
			a.change <- - d1 / d2
			a.change <- ifelse( abs( a.change ) > .3 , .3*sign(a.change) , a.change )              
			a.change <- a.change * est.aa
			a <- a + a.change
			a[ a < 0 ] <- .01			
#			cat( aa , " ") ; 
		if (progress){
				cat( paste( rep( "-" , prbar[aa]), collapse="") )
				flush.console()										
							
		if ( length(aG) < 2 ){ cat( paste( rep( "-" , 10 - length(aG) ), collapse="") ) }
						}
					}
		a1a <- max( abs( a - a0 ) )
		if (progress){ 
			cat("|     max. parm. change" , round( a1a , 5),"\n") 
				}
# cat("end a") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	

		#******************************************************************************
		# estimation of delta parameters
		delta0 <- delta
		dG <- setdiff( unique( est.delta ) , 0 )
		prbar <- seq( 1 , 10 , len = length(dG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		if (progress){ cat(" Estimation of delta: |") }
		if ( length(dG) == 0 ){ cat( paste(rep("-",10),collapse="") ) }
		if ( length(dG) > 0 ){
		# identify different a parameter groups
			est.cc <- 1 # * ( est.delta == cc )
			rescop <- .ll.rasch.copula320( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
								CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , 
								delta , wgt.theta , I , 
								bdat2.li , bdat2.li.resp  , pattern, GG , copula.type, Ncat.ld ,
							Qmatrix=Qmatrix		)
#			ll0 <- rescop$ll 
			rest1 <- .update.ll.rasch.copula321( theta.k , b, alpha1 , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta + h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type ,
							Qmatrix = Qmatrix )
			rest2 <- .update.ll.rasch.copula321( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta - h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type,
							Qmatrix=Qmatrix )
			ll0 <- ll1 <- ll2 <- rep(0,length(dG))
			rjkCC <- rescop$rjkCC
			for (cc in 1:CC){
					ll0[cc] <- sum( rjkCC[[cc]] * log( rescop$pjk.theta.kCC[[cc]] ) )
					ll1[cc] <- sum( rjkCC[[cc]] * log( rest1$pjk.theta.kCC[[cc]] ) )
					ll2[cc] <- sum( rjkCC[[cc]] * log( rest2$pjk.theta.kCC[[cc]] ) )				
								}
			a1 <- aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.delta) , sum , na.rm=T)				
			ll0 <- a1[,2]
			ll1 <- a1[,3]
			ll2 <- a1[,4]		
			delta.change <- nr.numdiff( ll0=ll0 , ll1=ll1 , ll2=ll2 , h=h )			

			ct <- sapply( dG , FUN = function(dd){
					( copula.type[ est.delta == dd ] )[1] } )
			maxstep <- ifelse( copula.type=="bound.mixt" , .2 , .9 )
			# hstep <- maxstep^( log( 2+iter))
			hstep_delta <- hstep <- hstep_delta * ( 1 / increment.factor )			
			
			delta.change <- ifelse( abs( delta.change ) > hstep , 
							hstep*sign(delta.change) , delta.change )              														
			delta.change <- delta.change[ match( est.delta , a1[,1] ) ]
			delta <- delta + delta.change
			delta[ delta <= 0 ] <- 2*numdiff.parm		
			delta <- ifelse( copula.type == "bound.mixt" & ( delta > 1 ) ,
						1 - 2 * numdiff.parm , delta )
			if (progress){ 
				cat( paste( rep( "-" , 10 ), collapse="") ) 
				flush.console()
						}
							}
        a1d <- max( abs( delta - delta0 ) )
		if (progress){ cat("|     max. parm. change" , round( a1d , 5),"\n") }
# cat("end delta") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
		#******************************************************************************
		# estimation of alpha parameters
		alpha10 <- alpha1
		alpha20 <- alpha2
		prbar <- 5
		if (progress){ cat(" Estimation of alpha: |")		 }
		# alpha1
		if (est.alpha){
			rescop <- .ll.rasch.copula320( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld ,
							Qmatrix=Qmatrix )
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula320( theta.k , b, alpha1 + h , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type,
							Qmatrix=Qmatrix )$ll
			ll2 <- .update.ll.rasch.copula320( theta.k , b , alpha1 - h , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG, copula.type,
							Qmatrix=Qmatrix )$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k1 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha1 <- alpha1 + alpha.change
			}
			if (progress){ 
				cat( paste( rep( "-" , prbar), collapse="") )
				flush.console()		
						}
		# alpha2
		if (est.alpha){
			rescop <- .ll.rasch.copula320( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG , copula.type, Ncat.ld,
							Qmatrix=Qmatrix )
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula320( theta.k , b, alpha1 , alpha2+h , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG , copula.type,
							Qmatrix=Qmatrix )$ll
			ll2 <- .update.ll.rasch.copula320( theta.k , b , alpha1 , alpha2 -h, a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG , copula.type,
							Qmatrix=Qmatrix )$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k2 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha2 <- alpha2 + alpha.change
			}			
			if (progress){
				cat( paste( rep( "-" , prbar), collapse="") )
				flush.console()		
						}
        a1k <- max( abs( c( alpha1 - alpha10, alpha2 - alpha20 )) )
		if (progress){ cat("|     max. parm. change" , round( a1k , 5),"\n")	}
		#******************************************************************************
		# estimation of mu parameters
		
		a1m <- 0
		if (G>1){
		mu0 <- mu
		# identify different a parameter groups
		muG <- 1:(GG-1)
		prbar <- seq( 1 , 10 , len = length(muG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		if (progress){ cat(" Estimation of mu:    |") }
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type , Ncat.ld)
			ll0 <- rescop$ll
			# mu + h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg] + h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld)
			ll1 <- rescop$ll
			# mu - h		
			w1 <- wgt.theta
			w2 <- stats::dnorm( theta.k , mean = mu[gg] - h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld )
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
			if (progress){
				cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
				flush.console()
						}
							}
		if ( length(muG) < 2 ){ cat( paste( rep( "-" , 10 - length(muG) ), collapse="") ) }
		a1m <- max( abs( mu - mu0 ) )
		if (progress){ cat("|     max. parm. change" , round( a1m , 5),"\n") }
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
		if (progress){ cat(" Estimation of sigma: |") }
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG,copula.type,Ncat.ld )
			ll0 <- rescop$ll
			# sigma + h		
			w1 <- wgt.theta
			w2 <- stats::dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] +h)
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type,Ncat.ld)
			ll1 <- rescop$ll
			# sigma - h		
			w1 <- wgt.theta
			w2 <- stats::dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg]-h )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type,Ncat.ld)
			ll2 <- rescop$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			sigma.change <- - d1 / d2
			sigma.change <- ifelse( abs( sigma.change ) > .3 , .3*sign(sigma.change) , sigma.change )      		
			sigma.change <- sigma.change * ( ( 1:G ) == gg )
			sigma <- sigma + sigma.change
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] )
			wgt.theta[,gg] <- w2 / sum(w2)
#			cat( aa , " ") ; 
		if (progress){
				cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
				flush.console()
						}
							}
		if ( length(sigmaG) < 2 ){ cat( paste( rep( "-" , 10 - length(sigmaG) ), collapse="") ) }
		a1s <- max( abs( sigma - sigma0 ) )
		if (progress){ cat("|     max. parm. change" , round( a1s , 5),"\n") }
				}			# end sigma
		######################################################################		
# cat("end rest") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1					
		# estimate covariance matrix		

		res0 <- rasch.copula3.covariance( f.qk.yi=rescop$post , Sigma=sigma , theta.k , N ,
					mu.fixed= mu.fixed , variance.fixed = NULL, D=D , est.corr=TRUE , irtmodel=NULL ,
					freqwgt = pattern$freqwgt)
		mu <- res0$mu
		sigma <- res0$Sigma
		sigma0 <- res0$Sigma0		
		wgt.theta <- res0$pi.k[,1]		
# cat("end covariance") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
		

#cat("other pars") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1						       				
		iter <- iter + 1 
#		thetawidth <- diff( theta.k )[1]		
		#**********************
		# deviance rasch.mml
        # ll[gg] <- sum( dat1[,2] * log( rowSums( f.yi.qk[,] * 
		#			outer( rep(1,nrow(f.yi.qk[,])) , pi.k ) ) ) )

		M2 <- outer( rep(1,nrow(pattern)), wgt.theta )	

		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm * M2 ) ) )
		#**********************
		# deviance tam
		# deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )		
#		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm) * thetawidth ) )
#		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm * M2 ) ))
        dev.change <- abs( ( dev - dev0)/ dev0 )
		absdev.change <- abs( dev- dev0 )
        par.change <- max( a1a , a1b , a1d , a1k , a1m , a1s)
		if (progress){
			cat( "Deviance = "  ,   round( dev , 5 ) , 
					" | Deviance change = " , - round( dev- dev0 , 4 ) ,
					"| max. parm. change = " ,  round( par.change , 6 ) ,  " \n"   )  
			if ( ( dev > dev0 ) & ( iter > 4 ) ){ cat("   Deviance has increased! Convergence Problems?\n") }
	#cat("rest") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1						       						
			flush.console()
					}				
			}
	# end MML iterations
	#**********************************************************************************
	# Standard error estimation (This is a TO DO!)
	iterend <- iter

	################################################
	# evaluation of posterior distribution
	# use functions from smirt
	
	post <- rescop$post
	EAP.Rel <- rep(0,D)
	names(EAP.Rel) <- paste0("Dim",1:D)
	for (dd in 1:D){
		#	dd <- 1
		M3 <- outer( rep(1,nrow(post)) , theta.k[,dd] ) 
		pattern[ , paste0("EAP.Dim",dd)] <- rowSums( M3 * post )				
		pattern[ , paste0("PostVar.Dim",dd)] <- 
			rowSums( M3^2 * post ) -  pattern[ , paste0("EAP.Dim",dd)]^2 				
		M.EAP <- weighted.mean( pattern[ , paste0("EAP.Dim",dd)] , pattern$freqwgt )
		Var.EAP <- sum(( ( pattern[ , paste0("EAP.Dim",dd)] - M.EAP )^2 * 
		                 pattern$freqwgt )) / ( sum( pattern$freqwgt ) )
		MVar.EAP <- weighted.mean( pattern[ , paste0("PostVar.Dim",dd)] , pattern$freqwgt )
		EAP.Rel[dd] <- Var.EAP / ( Var.EAP + MVar.EAP )
					}
	
	#********	

	
	#********************************************************
	# information criteria
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat00) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
		bG <- setdiff( unique( est.b ) , 0 )
		aG <- setdiff( unique( est.a ) , 0 )		
		dG <- setdiff( unique( est.delta ) , 0 )				
        ic[[ "np" ]] <- length(bG) + length(aG) + length(dG) + 2*est.alpha + D*(D-1)/2
        # AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np	
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )			
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
	if (progress){ 
		cat("--------------------------------------------------------------------------------\n")
					}
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
	if (progress){
		cat("Parameter summary\n")
		.pr( item , digits=3 )		# print item statistics
		# dependency parameter
		cat("\nDependency parameters\n")
				}
	summary.delta <- data.frame( "cluster" = 1:CC , "delta" = delta , 
				"est.delta" = est.delta , "copula.type" = copula.type	 )
	summary.delta$items <- sapply( 1:CC , FUN = function(cc){ 
		paste( colnames(dat)[ itemcluster == cc ] , collapse="-" )
				} )
			s2 <- Sys.time()				
	if (progress){ 				
		.pr(summary.delta , digits = 3)
		cat(paste("\nEAP Reliability:\n "))
		print( round( EAP.Rel,3))
		cat("\n")
		cat("Generalized logistic link function\n")
		cat("alpha1=",round(alpha1,3)," alpha2=" , round(alpha2,3) , " \n\n")	
			# computational time
		   cat("-----------------------------------------------------------------\n")
		   cat("Start:" , paste( s1) , "\n")
		   cat("End:" , paste(s2) , "\n")
		   cat("Difference:" , print(s2 -s1), "\n")
		   cat("-----------------------------------------------------------------\n")
				}
	datalist <- list( pattern.in.data = pattern.in.data , dat0 = dat0 ,
					dat2 = dat2 , dat2.resp = dat2.resp , dat2.li = dat2.li , 
					dat2.ld = dat2.ld , dat2.li.resp=dat2.li.resp , 
					dat2.ld.resp = dat2.ld.resp , dp.ld = dp.ld , CC = CC ,
					bdat2.li = bdat2.li , bdat2.li.resp = bdat2.li.resp , 
					itemcluster0 = itemcluster0 , dat3.ld = dat3.ld		
							)						
    # collect results
	v1 <- datalist$pattern.in.data	
	patternindex <- match( v1 , pattern$pattern )	
    names(mu) <- rownames(sigma) <- colnames(sigma) <- paste0("Dim",1:D)
	# person parameters
	person <- pattern[ patternindex, ] 
	person$freqwgt <- NULL		
	
	res <- list( "N.itemclusters" = CC , "item" = item , "iter" = iterend , "dev" = dev ,
					"delta" = delta , "b" = b , "a" = a , "mu" = mu , "sigma" = sigma , 
					"alpha1"=alpha1 , "alpha2"=alpha2 , "ic" = ic , "theta.k" = theta.k , 
					"pi.k"= wgt.theta , "deviance" = dev ,
					"pattern" = pattern	 , "person" = person ,
					"datalist" = datalist	, "EAP.Rel" = EAP.Rel	,
					"copula.type" = copula.type	, "summary.delta" = summary.delta ,
					"f.qk.yi" =(res.posterior$post)[ patternindex ,] , 
					"f.yi.qk" =(res.posterior$post.unnorm)[ patternindex ,] , 
					"Qmatrix" = Qmatrix , "D" = D , 					
					"s2" = s2 , "s1"=s1
								)	
	class(res) <- "rasch.copula3"
	return(res)
		}
#----------------------------------------------------------------------------------



