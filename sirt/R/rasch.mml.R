 

 
#------------------------------------------------------------------------
# --------------------------------------------------------------------#
# Semiparametric Maximum Likelihood Estimation in the Raschtype
#  Model                                                                            #
# item discrimination and guessing parameter can be fixed                                                                                    #
##NS export(rasch.mml)
rasch.mml <- function( dat , theta.k = seq(-4,4,len=20) , group = NULL , weights = NULL ,  
				constraints = NULL , 
                        glob.conv = 10^(-5) , parm.conv = 10^(-4) , mitermax = 1 , WLE = T , 
                        mmliter = 600 , progress = TRUE ,  
                        fixed.a = rep(1,ncol(dat)) , 
                        fixed.c = rep(0,ncol(dat)) , 
                        fixed.d = rep(1,ncol(dat)) , 
						fixed.K = rep(3,ncol(dat)) , 
						b.init = NULL , 
						est.a = rep(1,ncol(dat)) , est.b = NULL , 
						est.c = NULL , est.d = NULL , 
						max.c = 1 , min.d = 0 , est.K = NULL , max.K = 20 , 
                        pid = 1:( nrow(dat) ) , trait.weights = NULL ,
                        center.trait = TRUE , nplausible = 5 , alpha1 = 0 , alpha2 = 0 ,
                        est.alpha = FALSE , equal.alpha = FALSE , 
                        designmatrix = NULL , alpha.conv = parm.conv , 
						numdiff.parm = 0.00001 , numdiff.alpha.parm= numdiff.parm , 
						normal.trait = TRUE , ramsay.qm = FALSE ,  ... ){
    #******************************************************************************************##
    # INPUT:                                                                                ***##
    # dat           ... data frame with item responses
    # theta.k       ... grid of theta values where the trait density is evaluated
    # group         ... vector of group entries (numbered from 1 to G)
    # weights       ... sample weights or absolute frequency of item response patterns
    # glob.conv     ... global convergence criterion
    # conv1         ... convergence of parameters
    # mitermax      ... maximum number of iterations within the M step
    # WLE           ... should WLEs be estimated?
    # mmliter       ... maximum number of iterations during Maximum Likelihood Estimation
    # progress      ... should the estimation progress being displayed?
    # constraints   ... matrix with two columns: first column is item label (name or number), 
	#						second column is fixed item difficulty value (b parameter)
    # fixed.a       ... vector of fixed item discriminations
    # fixed.c       ... vector of fixed guessing parameters
    # b.init        ... initial estimates of item difficulties
	# est.a			... estimation of discrimination parameter
	# est.c, est.d	... estimated groups of parameters for lower and upper asymptote
	# max.c, min.d  ... maximal c and minimal d parameters to be estimated
    # pid           ... labels for subject ID's (-> pid ... person IDs)
    # trait.weights ... vector of fixed trait distribution
    # center.trait  ... set the mean of trait distribution equal to zero? (default = TRUE)
    # nplausible    ... number of plausible values to be drawn for item fit estimation
    # est.alpha     ... should alpha parameters be estimated?
    # equal.alpha   ... should equal alpha's be estimated?
    # designmatrix  ... Q matrix of item parameter restrictions    
    # numdiff.parm  ... step parameter for numerical differentiation
    # normal.trait  ... normal distribution assumption of the trait
	# ramsay.qm		... estimate quotient model of ramsay?
    #******************************************************************************************##
    #****
    # specifications
    conv1 <- parm.conv
#	ramsay.qm <- FALSE
    if ( ! ramsay.qm) { 
		raschtype <- TRUE  
			}
	if (ramsay.qm ){
		raschtype <- FALSE
		# no alpha, a, c or d parameters can be estimated
		est.alpha <- FALSE	
		est.a <- est.c <- est.d <- NULL	
		pow.qm <- 1	# This parameter is ignored in analyses
					}
    # computation time
    s1 <- Sys.time()
    if (est.alpha){ 
            if (is.null(alpha1) ){ alpha1 <- 0 }
            if (is.null(alpha2) ){ alpha2 <- 0 }
                }
    #**************************************************************************************
    # some data checks
	ag1 <- NULL
	if( max( colMeans( is.na( dat ) ) ) == 1 ){
		stop("Remove items which have no observations!")
					}
    if ( ! is.null(group) ){ 
            t1 <- table(sort(group) )
			group.orig <- group
			group <- match( group.orig , sort(unique( group.orig)) )	
			ag1 <- stats::aggregate( group , list( group.orig) , mean )
			colnames(ag1) <- c("group" , "groupindex" )
#            x1 <- seq( 1 , length(t1) ) 
#            x2 <- as.numeric(sort(names(t1)) )
#            if ( sum( abs( x1-x2) ) > 10^(-10) ){
#                stop( "Group identifier must be recoded to 1, ..., G" )
#                                    }
                            }
    # center trait: if there exists constraints, then do not center
	
	if ( is.null( colnames(dat) ) ){ 
			colnames(dat) <- paste( "I" , 1:ncol(dat) , sep="")
						}
	  if ( ! is.null( constraints ) ){ 
            center.trait <- F 
          if( ! is.numeric( constraints[,1] ) ){
            constraints[,1] <- match( paste(constraints[,1]) , colnames(dat) )
                                       }           
             constraints <- na.omit(constraints)
             constraints <- constraints[ constraints[,1] <= ncol(dat) , ]
                                }
    if ( ! is.null( designmatrix) ){
            if ( ncol(dat) != nrow(designmatrix) ){ 
                    stop( "Row dimension of designmatrix should be equal to number of items")
                                }
                        }
	# est.b parameters
	if (! is.null(est.b) ){
		bG <- unique( est.b ) 
		designmatrix <- matrix( 0 , ncol(dat) , length(bG) )
		for (bb in bG){
			# bb <- bG[1]
			designmatrix[ which( est.b == bb ) , bb ] <- 1
						}
				}
	# set starting values for estimated c and d parameters
	if ( sum(est.c) > 0 ){	fixed.c[ est.c > 0 ] <- .10 }
	if ( sum(est.d) > 0 ){	fixed.d[ est.d > 0 ] <- .95 }

	
    #****************************************************************************************
     WLE <- FALSE 
     pure.rasch <- -9   # this parameter is only included for historical reasons of this program.
    # specify weights
    if ( is.null(weights) ){ weights <- rep( 1 , nrow(dat) ) }
    # display
    if ( progress & ( ramsay.qm ) ){
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
        cat("Quotient Model (Ramsay, 1989) \n") 
        if (normal.trait){ cat("Normal trait distribution \n") } else { cat("Nonparametric trait distribution \n") }
#		if (ramsay.qm){ cat("Log Normal Distribution of Theta with Power of" , pow.qm , "\n") }
        cat("---------------------------------------------------------------------------------------------------------- \n")
        flush.console()
      }
    if ( progress & (raschtype) ){
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
        if ( est.alpha ){ 
            cat(paste( "Raschtype Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
                    } else {
            cat(paste( "Raschtype Model with generalized logistic link function: alpha1=",alpha1 , " , alpha2=", alpha2 , " \n") )
                        }
		if ( sum(est.c) > 0){ cat(paste( "Estimated guessing parameter groups \n") )}  ## estimated guessing parameters
		if ( sum(est.d) > 0){ cat(paste( "Estimated slipping parameter groups \n") )}  ## estimated slipping parameters
        cat("---------------------------------------------------------------------------------------------------------- \n")
        flush.console()
      }

    # revise guessing parameter (if necessary)
    if ( !is.null(fixed.c) ){
        # calculate itemmeans
        itemmean <- colMeans( dat ,  na.rm = T )
        if ( any( itemmean < fixed.c) ){  
                cat ( "revise fixed guessing estimates\n")
                fixed.c[ itemmean < fixed.c] <- 0
                }
            }
        # data preparations 
        dp <- .data.prep( dat , weights = weights , use.freqpatt = ( is.null(group) ) )
        dat1 <- dp$dat1
        dat2 <- dp$dat2
        dat2.resp <- dp$dat2.resp
        freq.patt <- dp$freq.patt
        n <- dp$n
        I <- dp$I
        # probability weights at theta.k
        pi.k <- stats::dnorm( theta.k ) 
        pi.k <- pi.k / sum( pi.k )
        # group calculations
        if ( !is.null( group )){ 
            G <- length( unique( group ) )
            pi.k0 <- pi.k
            pi.k <- matrix( 0 , nrow=length(pi.k0) , ncol=G)
            for (gg in 1:G){
                pi.k[,gg] <- pi.k0
                              }
                    }
        # initial estimates for item difficulties
        if ( is.null(b.init) ){   
			b <- - stats::qlogis( colMeans( dat , na.rm=T ) )
				if ( FALSE ){ 			
#				if ( ramsay.qm ){ 			
						b <-   - log( ( fixed.K * colSums( dat , na.rm=TRUE ) ) / 
									( colSums( 1 - dat , na.rm=TRUE ) ) ) 
								}
						} else { 
					b <- b.init 
					}
        # initial iteration index
        iter <- 0
        par.change <- dev.change <- 3
        dev <- 99
        maxalphachange <- 1
		  # display
		  disp <- "...........................................................\n"		
        #--------------------------------#
        # MML Iteration Algorithm        #
        while ( ( dev.change > glob.conv | par.change > conv1 | maxalphachange > alpha.conv ) & iter < mmliter ){
		if (progress){ 
		  cat(disp)	
		  cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )
          flush.console()
		}		
				b0 <- b
                dev0 <- dev    		
                # perform E Step
                if ( ! raschtype ){ 
					e1 <- .e.step.ramsay( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
									fixed.K , group , pow.qm = pow.qm  )
                                } else { 
                    e1 <- .e.step.raschtype( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
                                    fixed.a , fixed.c , fixed.d ,  alpha1 , alpha2 , group )        
                                }
                n.k <- e1$n.k
                n.jk <- e1$n.jk
                r.jk <- e1$r.jk
                pjk <- e1$pjk
                f.qk.yi <- e1$f.qk.yi
                dev <- -2*e1$ll		
                # perform M Step
                 if ( ! raschtype ){ 
                    m1 <- .m.step.ramsay( theta.k , b , n.k , n , n.jk , r.jk , I , conv1 , constraints , 
							mitermax , pure.rasch ,  trait.weights , fixed.K , 
							designmatrix = designmatrix , group = group ,                                            numdiff.parm=numdiff.parm , pow.qm = pow.qm )
											} else {
                    m1 <- .m.step.raschtype( theta.k , b , n.k , n , n.jk , r.jk , 
							I , conv1 , constraints , mitermax , pure.rasch ,                                             trait.weights , fixed.a , fixed.c , fixed.d ,  alpha1 , 
							alpha2 , designmatrix = designmatrix ,
							group = group ,
                            numdiff.parm=numdiff.parm)
                                }
 								
                b <- m1$b
                # distribution
                G <- m1$G
                pi.k <- m1$pi.k
                # latent ability distribution
                if (normal.trait){ 
                    sd.trait <- mean.trait <- rep(0,G)
                    for (gg in 1:G){ 
                        mean.trait[gg] <- stats::weighted.mean( theta.k , pi.k[,gg] )
                        sd.trait[gg] <- sqrt( stats::weighted.mean( ( theta.k - mean.trait[gg] )^2 , pi.k[,gg] ) ) 
                        if (center.trait){ mean.trait[1] <- 0 } 
						if ( ! is.null(est.a) ){ sd.trait[1] <- 1 }	
                        pi.k[,gg] <- stats::dnorm( theta.k , mean = mean.trait[gg] , sd = sd.trait[gg] )
                        pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
							}
					# sigma <- sd.trait
						}
        ##############################
        # estimation of alpha, c and d parameters
        alpha.change <- 0
        maxalphachange <- 0
		a1a <- a1b <- 0
		a1K <- a1c <- 0
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of a parameters
		if ( sum( est.a ) > 0 ){
                h <- numdiff.parm 
				fixed.a0 <- fixed.a
				# identify different a parameter groups
				aG <- setdiff( unique( est.a ) , 0 )
				for (aa in aG){
				# cc <- cG[1]	
					est.aa <- 1 * (est.a == aa )
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c  )
					fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d )            
					#****
					# alpha1
					pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a)
						pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a + h*est.aa)
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c )				
						pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a - h*est.aa)
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c )				
						pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
					# first order derivative
					# f(x+h) - f(x-h) = 2* f'(x) * h
					ll0 <- ll1 <- ll2 <- rep(0,G)
				# was ist hier das G? => G ist hier Anzahl der Gruppen
					for (gg in 1:G){ 
						ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
								}  
					ll0 <- sum(ll0)
					ll1 <- sum(ll1)
					ll2 <- sum(ll2)								
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
#					se.a1[aa] <- sqrt( -1 / d2  )
					# change in item difficulty
					a.change <- - d1 / d2
					a.change <- ifelse( abs( a.change ) > .1 , .1*sign(a.change) , a.change )              
					a.change <- a.change * est.aa
					fixed.a <- fixed.a + a.change
					fixed.a[ fixed.a < 0 ] <- 0
					# fixed.a[ fixed.a > max.c ] <- max.c
								}
                a1a <- max( abs( fixed.a - fixed.a0 ) )
						}

						
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of c parameters
		if ( sum( est.c ) > 0 ){
                h <- numdiff.parm 
				fixed.c0 <- fixed.c
				# identify different c parameter groups
				cG <- setdiff( unique( est.c ) , 0 )
				for (cc in cG){
				# cc <- cG[1]	
					est.cc <- 1 * (est.c == cc )
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c  )
					fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d )            
					#****
					# alpha1
					pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a)
						pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a)
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c + h*est.cc )				
						pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a)
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c - h*est.cc )				
						pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
					# first order derivative
					# f(x+h) - f(x-h) = 2* f'(x) * h
					ll0 <- ll1 <- ll2 <- rep(0,G)
					for (gg in 1:G){ 
						ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
								}     
					ll0 <- sum(ll0)
					ll1 <- sum(ll1)
					ll2 <- sum(ll2)								
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
					# change in item difficulty
					c.change <- - d1 / d2
					c.change <- ifelse( abs( c.change ) > .1 , .1*sign(c.change) , c.change )              
					c.change <- c.change * est.cc
					fixed.c <- fixed.c + c.change
					fixed.c[ fixed.c < 2*numdiff.parm ] <- 2*numdiff.parm
					fixed.c[ fixed.c > max.c ] <- max.c
								}
                a1b <- max( abs( fixed.c - fixed.c0 ) )
						}
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of d parameters
		if ( sum( est.d ) > 0 ){
                h <- numdiff.parm 
				fixed.d0 <- fixed.d
				# identify different c parameter groups
				dG <- setdiff( unique( est.d ) , 0 )
				for (dd in dG){
				# cc <- cG[1]	
					est.dd <- 1 * (est.d == dd )
					fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c  )
					fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d )            
					#****
					# alpha1
					pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a)
						pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a)
					fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d + h*est.dd)     					
						pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a)
					fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d - h*est.dd)     					
						pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
					# first order derivative
					# f(x+h) - f(x-h) = 2* f'(x) * h
					ll0 <- ll1 <- ll2 <- rep(0,G)
					for (gg in 1:G){ 
						ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
								}  
					ll0 <- sum(ll0)
					ll1 <- sum(ll1)
					ll2 <- sum(ll2)								
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
					# change in item difficulty
					d.change <- - d1 / d2
					d.change <- ifelse( abs( d.change ) > .1 , .1*sign(d.change) , d.change )              
					d.change <- d.change * est.dd
					fixed.d <- fixed.d + d.change
					fixed.d[ fixed.d > 1 - 2 * numdiff.parm ] <- 1 - 2*numdiff.parm
					fixed.d[ fixed.d < min.d ] <- min.d
								}
                a1c <- max( abs( fixed.d - fixed.d0 ) )
						}
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of d parameters
		if ( sum( est.K ) > 0 & ramsay.qm ){
                h <- numdiff.parm 
				fixed.K0 <- fixed.K
				# identify different c parameter groups
				kG <- setdiff( unique( est.K ) , 0 )
				for (kk in kG){
				# cc <- cG[1]	
					est.kk <- 1 * (est.K == kk )
					#.prob.ramsay <- function( theta , b ,  fixed.K = 3 + 0*b)			
					#****
					# kk
					pjk <- .prob.ramsay( theta.k , b , fixed.K , pow = pow.qm )
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.ramsay( theta.k , b , fixed.K+h*est.kk , pow = pow.qm  )
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.ramsay( theta.k , b , fixed.K-h*est.kk , pow = pow.qm  )
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
					# first order derivative
					# f(x+h) - f(x-h) = 2* f'(x) * h
					ll0 <- ll1 <- ll2 <- rep(0,G)
					for (gg in 1:G){ 
						ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
								}  
					ll0 <- sum(ll0)
					ll1 <- sum(ll1)
					ll2 <- sum(ll2)								
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
					# change in item difficulty
					K.change <- - d1 / d2
					K.change <- ifelse( abs( K.change ) > 1 , 1*sign(K.change) , K.change )              
					K.change <- K.change * est.kk
					fixed.K <- fixed.K + K.change
					fixed.K[ fixed.K < 1 ] <- 1 
					fixed.K[ fixed.K > max.K ] <- max.K
#					fixed.K[ fixed.d < min.d ] <- min.d
								}
				# convergence is indicated in metric guess.K = 1 / ( fixed.K + 1 )
                a1K <- max( abs( 1/(1+fixed.K) - 1/(1+fixed.K0) ) )
						}						
		#***************************
		# estimation of alpha
        if ( est.alpha ){ 
                alpha1.old <- alpha1
                h <- numdiff.alpha.parm 
#			if (iter > 30 ){ 
#					h <- .01 
#					print(h)
#						}
                fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c )
                fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d )            
                #****
                # alpha1
                pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a)
                    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
                pjk <- ( pjk + .000000005 ) / 1.00000001 
                pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
                pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1+h , alpha2 , fixed.a)
                    pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
                pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
                pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
                pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1-h , alpha2 , fixed.a)
                    pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
                pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
                pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
                # first order derivative
                # f(x+h) - f(x-h) = 2* f'(x) * h
                ll0 <- ll1 <- ll2 <- rep(0,G)
                for (gg in 1:G){ 
                    ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
                    ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
                    ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
                            }                
                ll0a1 <- ll0 <- sum(ll0)
                ll1a1 <- ll1 <- sum(ll1)
                ll2a1 <- ll2 <- sum(ll2)
                d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
#				d1 <- ( ll1 - ll0 ) / h
                # second order derivative
                # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
                d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
                # change in item difficulty
                alpha.change <- - d1 / d2
                alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
                alpha1 <- alpha1 + alpha.change
                a1 <- abs(alpha.change )
                #****
                # alpha2
                pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a)
                    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
                pjk <- ( pjk + .000000005 ) / 1.00000001 
                pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
                pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 +h , fixed.a )
                    pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
                pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
                pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
                pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 -h , fixed.a)
                    pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
                pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
                pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
                # first order derivative
                # f(x+h) - f(x-h) = 2* f'(x) * h
                ll0 <- ll1 <- ll2 <- rep(0,G)
                for (gg in 1:G){ 
                    ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
                    ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
                    ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
                            }                
                ll0 <- sum(ll0)
                ll1 <- sum(ll1)
                ll2 <- sum(ll2)
                d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
#				d1 <- ( ll1 - ll0 ) / h

                # second order derivative
                # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
                d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
                # change in item difficulty
                alpha.change <- - d1 / d2
                alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )
                alpha2 <- alpha2 + alpha.change
                a2 <- abs(alpha.change)
                maxalphachange <- max(a1, a2)
                if (equal.alpha){
                        ll0 <- ll0a1 + ll0 
                        ll1 <- ll1a1 + ll1                 
                        ll2 <- ll2a1 + ll2                 
                        d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
                        d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2                
                        alpha.change <- - d1 / d2
                        alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )
                        alpha2 <- alpha1 <- alpha1.old + alpha.change
                        a2 <- abs(alpha.change)
                        maxalphachange <- max(a2)
                            }
						}
            ############################
                # possibly incorrect deviance calculated at the M step
                #                dev <- -2*m1$ll
                # iteration index
                dev.change <- abs( ( dev - dev0)/ dev0 )
                par.change <- max( c( abs(b - b0 ) , abs(alpha.change ) , a1a , a1b , a1c , a1K ) )
                # display convergence
                if ( progress  ){   
						cat( paste( "   Deviance = "  , 
                             round( dev , 4 ) , 
							 if (iter > 0 ){ " | Deviance change = " } else {""} ,
							if( iter>0){round( - dev + dev0 , 6 )} else { ""}	,						 
							 "\n    Max. parm. change = " , 
                             round( par.change , 6 ) ,  " \n"   )  )  
                                    if ( est.alpha ){           
                                        cat( paste( "    alpha1=" , round(alpha1,3) , " | alpha2=", round( alpha2,3) , 
                                                    " | max alpha change " , round( maxalphachange ,7 ) , "\n" , sep=""))
												}
                                    if ( sum(est.a) > 0  ){           
										apars <- stats::aggregate( fixed.a , list(est.a) , mean )
																apars <- apars[ apars[,1] > 0 , ]
																a.est <- apars[  , 2 ]
                                        cat( paste( "    Estimated discrimination parameter groups a = " , 
												paste( round(a.est ,3) , collapse=" " ) , "\n" , sep=""))
												}
									if ( sum(est.c) > 0  ){           
										cpars <- stats::aggregate( fixed.c , list(est.c) , mean )
																cpars <- cpars[ cpars[,1] > 0 , ]
																c.est <- cpars[  , 2 ]
                                        cat( paste( "    Estimated guessing parameter groups c = " , 
												paste( round(c.est ,3) , collapse=" " ) , "\n" , sep=""))
												}
                                    if ( sum(est.d) > 0  ){           
										dpars <- stats::aggregate( fixed.d , list(est.d) , mean )
																dpars <- dpars[ dpars[,1] > 0 , ]
																d.est <- dpars[  , 2 ]
                                        cat( paste( "    Estimated slipping parameter groups d = " , 
												paste( round(d.est ,3) , collapse=" " ) , "\n" , sep=""))
												}
                                    if ( sum(est.K) > 0  ){           
										dpars <- stats::aggregate( fixed.K , list(est.K) , mean )
																dpars <- dpars[ dpars[,1] > 0 , ]
																K.est <- dpars[  , 2 ]
                                        cat( paste( "    Estimated K parameter groups K = " , 
												paste( round(K.est ,3) , collapse=" " ) , "\n" , sep=""))
												}												
									flush.console() 
                                    }
                iter <- iter + 1
                     }
        #**********************************************
        # standard error for item parameter
		# 		THIS IS NOT YET DONE!
		# ...
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
        ic[[ "np" ]] <- G + ncol(dat) + ( G - 1 )
        # subtract fixed constraints
        if ( ! is.null( constraints) ){ ic$np <- ic$np - nrow(constraints) }
        # subtract constraints due to designmatrix
        if ( ! is.null( designmatrix ) ){ ic$np <- ic$np - ncol(dat) + ncol(designmatrix) }
        # alpha estimation
        ic$np <- ic$np + est.alpha * 2 - equal.alpha *1
		# guessing, slipping and discrimination parameter estimation
		if ( sum(est.c) > 0 ){ ic$np <- ic$np + length(cG) }
		if ( sum(est.d) > 0 ){ ic$np <- ic$np + length(dG) }		
		if ( sum(est.a) > 0 ){ ic$np <- ic$np + length(aG) }		
		if ( sum(est.K) > 0 ){ ic$np <- ic$np + length(kG) }		
        # AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
        # item statistics
        item <- data.frame( "item" = colnames(dat) , "N" = colSums( weights*(1 - is.na(dat)) ) , 
                            "p" = colSums( weights*dat   , na.rm=T) / colSums( weights*(1 - is.na(dat)) ) , 
                            "b" = b  )
		if (is.null(est.b)){ item$est.b = seq(1,I) } else { item$est.b <- est.b }
        # fixed parameters
        item$a <- fixed.a
		if ( ! is.null( est.a) ){  item$est.a <- est.a } else { item$est.a <- rep(0,I) }
		# include threshold
		item$thresh <- item$a*item$b		
		# guessing parameter
        item$c <- fixed.c
		if ( ! is.null( est.c) ){ item$est.c <- est.c } else { item$est.c <- rep(0,I) }
        item$d <- fixed.d
		if ( ! is.null( est.d) ){  item$est.d <- est.d } else { item$est.d <- rep(0,I) }
        if (m1$center){  if ( is.null(constraints) ){ #    item[I,4] <- NA 
                                                                    } 
                    else { item[ constraints[,1] ,4] <- NA } }
        rownames(item) <- colnames(dat)
        # latent ability distribution
        sd.trait <- mean.trait <- rep(0,G)
        for (gg in 1:G){ 
            mean.trait[gg] <- weighted.mean( theta.k , pi.k[,gg] ) 
            sd.trait[gg] <- sqrt( weighted.mean( ( theta.k - mean.trait[gg] )^2 , pi.k[,gg] ) ) 
                    }
        # center trait distribution
        if ( center.trait & G < 1 ){
                theta.k <- theta.k - mean.trait
                b <- b - mean.trait
                item$itemdiff <- b
                mean.trait <- 0
                        }      
        trait.distr <- data.frame( "theta.k" = theta.k , "pi.k" = pi.k )
        # item response pattern
        ability.est <- data.frame( dat1 , theta.k[ whichrowMaxs( f.qk.yi )$arg ] )
        colnames(ability.est) <- c("pattern" , "AbsFreq" , "mean" , "MAP" )
        ability.est$EAP <- rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , theta.k  )  )
        ability.est$SE.EAP <- sqrt( rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , theta.k^2  )  ) - ability.est$EAP^2 )
        # posterior distribution
        rownames(f.qk.yi) <- dat1[,1]
        # merging ability estimates
      #   
        if ( ! is.null(group)){  
                    ability.est2 <- cbind( freq.patt , ability.est[,-1] ) 
                            } else {
                ability.est2 <- merge( freq.patt , ability.est , 1 , 1 )
                            }
        ability.est2 <- ability.est2[ order(ability.est2$index) , -c(3:5) ]   
		# EAP reliability estimate
		reliability <- NULL
		reliability$eap.reliability <- 1 - mean(ability.est2$SE.EAP^2) / ( mean(ability.est2$SE.EAP^2) + var(ability.est2$EAP) )
        # include person ID
        ability.est2$pid <- pid
        # output fixed.a and fixed.c         
        if ( is.null(fixed.a ) & is.null(fixed.c) ){  fixed.a <- rep(1,I) ; fixed.c <- rep(0,I) }
        # include item discrimination
        item$emp.discrim <- round( item.discrim( dat ,  ability.est2$MAP ) , 3 )
		item$alpha1 <- alpha1
		item$alpha2 <- alpha2
		#---------------------------------------------------------
		# item summary Ramsay QM
		item2 <- NULL
		if ( ramsay.qm){
		if ( is.null(est.K) ){ est.K <- rep(0,I) }
			item2 <- data.frame( "item" = item$item  , "N" = item$N	 , "p" = item$p , 
							"K" = fixed.K , "est.K" = est.K , 
							"b" = exp(b) , "log_b" = b , "est.b" = item$est.b , 
							"guess.K" = 1/(fixed.K+1) , 
							"emp.discrim" = item$emp.discrim )
							}
        # result
        res <- list( "dat" = dat , "item" = item , "item2" = item2 , "trait.distr" = trait.distr , "mean.trait" = mean.trait , "sd.trait" = sd.trait , 
                    "deviance" = dev ,    "person" = ability.est2 , "pid" = pid , 
                    "ability.est.pattern" = ability.est , "f.qk.yi" =  f.qk.yi ,
                    "pure.rasch" = pure.rasch  , "fixed.a" = fixed.a , "fixed.c" = fixed.c , "G" = G ,
                    "alpha1"=alpha1 , "alpha2" = alpha2 , "iter" = iter , "reliability" = reliability ,
					"ramsay.qm" = ramsay.qm ) 
        class(res) <- "rasch.mml"
        # estimate plausible values
#        res$plausible.values <- plausible.values( res , nplausible = nplausible )	
# eliminate plausible value imputation 2012-02-17
        res$ic <- ic
		res$est.c <- est.c
		res$groupindex <- ag1
        # computation time
        s2 <- Sys.time()
		res$s1 <- s1
		res$s2 <- s2
		res$Rfcttype <- "rasch.mml"
        if (progress){ 
                cat("---------------------------------------------------------------------------------------------------------- \n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("---------------------------------------------------------------------------------------------------------- \n")
                    }       
        #..................
        return( res )
        }
#----------------------------------------------------------------------





