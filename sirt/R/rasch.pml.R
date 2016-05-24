

#####################################################
# Pairwise marginal likelihood (PML) estimation
##NS export(rasch.pml)
rasch.pml <- function( dat , est.b = seq( 1 , ncol(dat) ) , 
			est.a = rep( 0 , ncol(dat) ) , 
			est.sigma = TRUE ,
			itemcluster = NULL , 
#			Q=NULL ,		
#			zero.corrs= NULL , 
			weight = rep(1,nrow(dat)) ,  
			numdiff.parm=.001 , b.init = NULL , 
			a.init=NULL , sigma.init = NULL , 
			error.corr = 0*diag( 1 , ncol(dat) ) ,
            glob.conv=.001 , conv1 = .001 , pmliter = 100 ,
			progress = TRUE  ){
	##################################
	# multidimensional version does not work
	Q <- NULL ; combs <- zero.corrs <- NULL
	##################################
    # load libraries
    # extract information from data
	link <- "probit"
	s1 <- Sys.time()
	I <- ncol(dat)
	if( is.null( colnames(dat)) ){ colnames(dat) <- paste("I" , 1:I , sep="") }
	dat <- dat[ rowMeans( is.na( dat) ) < 1 , ]
	dat0 <- dat
	dat[ is.na(dat) ] <- 9
    N <- nrow(dat)
	    if ( progress  ){
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Pairwise Marginal Likelihood Estimation \n")
        cat(paste( "Raschtype Model with" , link , "link" ) , "\n") 
        cat("---------------------------------------------------------------------------------------------------------- \n")
        flush.console()
      }
    # extract frequencies from items
    p1 <- t(( dat == 1 )) %*% ( dat == 1 )
    itemfreq <- data.frame( "item" = 1:I , "p1" = diag(p1)  )
    # extract frequencies from item pairs
    itempairs <- as.data.frame( t( combn(1:I , 2 ) ) )
    colnames(itempairs) <- c("item1" , "item2")
    itempairs <- itempairs[ order( 1000 * itempairs[,2]  + itempairs[,1]  ) , ]
	# create weight matrix
	weight <- nrow(dat) * weight / sum(weight)
#	weight <- rep(1,nrow(dat))
	WM <- sqrt( outer( weight , rep(1,I) ) )
    # p11
    p11 <- t(( dat == 1 )*WM) %*% (( dat == 1 )*WM )
    itempairs$f11 <- p11[ upper.tri(p11) ]
    # p10
    p10 <- t(( dat == 1 )*WM) %*% (( dat == 0 )*WM )
    itempairs$f10 <- p10[ upper.tri(p10) ]
    # p01
    p10 <- t(( dat == 0 )*WM) %*% (( dat == 1 )*WM )
    itempairs$f01 <- p10[ upper.tri(p10) ]
    # p01
    p10 <- t(( dat == 0 )*WM) %*% (( dat == 0 )*WM )
    itempairs$f00 <- p10[ upper.tri(p10) ]
	
	# error correlations
	est.eps.corr <- itempairs$est.eps.corr <- error.corr[ upper.tri( error.corr ) ]
	est.corrs <- FALSE
	if ( any( itempairs$est.eps.corr != 0 ) ){ est.corrs <- TRUE }
	if ( link == "logit" ){ est.corrs <- FALSE }
	eps.corr <- 0.2 * ( itempairs$est.eps.corr != 0 )
	
	#*********
	# exclude some item pairs from calculation because they are
	# located in the same itemclusters (local dependence)
	# if error calculations are estimated, then no itemclusters 
	# can be selected
	if ( sum( error.corr) > 0 ){ itemcluster <- NULL }
	
	itemclusters <- unique( itemcluster[ itemcluster != 0 ] )
	IC <- length( itemclusters)
	for ( cc in itemclusters ){	
	#	cc <- 1
		icc <- which( itemcluster == cc )
		elim <- intersect( which( itempairs[ , "item1" ] %in% icc ) , which( itempairs[ , "item2" ] %in% icc ) )
		itempairs <- itempairs[ - elim , ]
				}
	if ( ! is.null( itemcluster ) ){
		eps.corr <- rep(0,nrow(itempairs ))
							}
    #*******
    # evaluate pairwise likelihood
	if ( is.null(b.init)){ 
		b <- b0 <- - stats::qnorm( colMeans( dat0 , na.rm=T) )
     if ( sum( est.b != 	seq( 1 , ncol(dat) ))>0 ){
		b <- 0*b
				}
			} #	if ( link == "logit"){
#		b <- b0 <- - qlogis( colMeans( dat0 , na.rm=T) )
#		p.ki <- c( .25220 , .58522 , .16257 )
#		S.ki <- c( .90793 , .57778 , .36403 )
#				}
	if ( sum(est.a) >0 ){ a <- rep(.5,I) } else { a <- rep(1,I) }
	a1s <- a1b <- a1a <- 0
	if ( ! is.null( b.init) ){ b <- b0 <- b.init }
	if ( ! is.null( a.init) ){ a <- a0 <- a.init }	
    if ( is.null(sigma.init)){ sigma <- 1 } else { sigma <- sigma.init }
	#******
				D <- 1 
#				}
    IP <- nrow(itempairs )
    itempairs$p1.item2 <- itempairs$p1.item1 <- rep(0,IP)
    itempairs$p11 <- rep(0,IP)
    iter <- 0
    #**********************************
    # BEGIN MARGINAL MAXIMUM LIKELIHOOD ESTIMATION
    dev <- 1 ; par.change <- dev.change <- 1000 
    while ( ( dev.change > glob.conv | par.change > conv1  ) & ( iter < pmliter )	){        
        cat( paste(rep("-" , 70), collapse="") , "\n")
        k1 <- floor( log10(iter+1) )
        x1 <- "        |" 
        x1 <- substring( x1 , k1+1 )
        s1c <- Sys.time()
        cat( paste( paste( "PML EM Iter." , iter + 1 ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
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
            # bb <- bG[1]
#            est.bb <- est.b * (est.b == bb )
			est.bb <- 1*(est.b == bb )
				if (bb == 1 ){ 
					respml0 <- .ll.rasch.pml.probit( b ,a ,  sigma , Q , 
							eps.corr ,itempairs  , IP , eps=10^(-14) )
					ll0 <- respml0$ll
								}                            
				 respml <- .update.ll.rasch.pml.probit( respml0 , b0 + h*est.bb , a , sigma , 
							Q , eps.corr ,itempairs  , IP , eps=10^(-14) )
				 ll1 <- respml$ll
				 respml <- .update.ll.rasch.pml.probit( respml0 , b0 - h*est.bb , a , sigma , 
								Q , eps.corr ,itempairs  , IP , eps=10^(-14) )
				 ll2 <- respml$ll                           
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2               
            if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
            b.change <- - d1 / d2
            b.change <- ifelse( abs( b.change ) > .5 , .5*sign(b.change) , b.change )             
            b.change <- b.change * est.bb
            b <- b + b.change
#           cat( bb , " ") ; 
            cat( paste( rep( "-" , prbar[match(bb,bG)] ), collapse="") )
            flush.console()             
                    }
        a1b <- max( abs( b - b0 ) )
        cat("|     max. parm. change" , round( a1b , 5),"\n")
		

        #************************************
        # estimation of a parameters
        a0 <- a
        # identify different a parameter groups
        aG <- setdiff( unique( est.a ) , 0 )
        prbar <- seq( 1 , 10 , len = length(aG) )
        prbar <- floor( prbar )
        prbar <- c( prbar[1] , diff(prbar) )
		if ( sum(est.a) > 0 ){
			cat(" Estimation of a:     |")
			for (aa in aG){
				# bb <- bG[1]
				est.aa <- 1*(est.a == aa )      
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma , 
					Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll0 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a+h*est.aa , sigma  , 
						Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll1 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a-h*est.aa , sigma , 
						Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll2 <- respml0$ll							 
				d1 <- ( ll1 - ll2  ) / ( 2 * h )    
				# second order derivative
				# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
				d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2               
				if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
				a.change <- - d1 / d2
				a.change <- ifelse( abs( a.change ) > .2 , .2*sign(a.change) , a.change )             
				a.change <- a.change * est.aa
				a <- a + a.change
				a[ a < .01 ] <- .01
				cat( paste( rep( "-" , prbar[match(aa,aG)] ), collapse="") )
				flush.console()             
						}
			a1a <- max( abs( a - a0 ) )
			cat("|     max. parm. change" , round( a1a , 5),"\n")		
					}		
		######################################
		# estimation sigma
		sigma0 <- sigma
		cat(" Estimation of sigma: |")	
		if (est.sigma & is.null(Q) ){
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma , 
								Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll0 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma + h , 
								Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll1 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma - h, 
							Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll2 <- respml0$ll		
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k2 <- alpha.change <- ifelse( abs( alpha.change ) > .2 , .2*sign(alpha.change) , alpha.change )              
			sigma <- sigma + alpha.change
								}
			prbar <- 10 
			flush.console()	
		if (est.sigma & (!is.null(Q)) ){
			for ( zz in seq(1,nrow(combs)) ){
				ii <- combs[zz,1]
				jj <- combs[zz,2]
				sigma.hh <- 0*sigma
				sigma.hh[ii,jj] <- sigma.hh[jj,ii] <- 1
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma , 
								Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll0 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma + h*sigma.hh , 
								Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll1 <- respml0$ll
				  respml0 <- .ll.rasch.pml.probit( b , a , sigma - h*sigma.hh, 
							Q , eps.corr , itempairs  , IP , eps=10^(-14) )
				  ll2 <- respml0$ll		
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
			alpha.change <- - d1 / d2
			a1k2 <- alpha.change <- ifelse( abs( alpha.change ) > .12 , .12*sign(alpha.change) , alpha.change ) 
			alpha.change <- alpha.change*sigma.hh
		
			sigma <- sigma + alpha.change
									}
		diag(sigma)[ diag(sigma) < .001 ] <- .0001
								}				
		cat( paste( rep( "-" , prbar), collapse="") )
        a1s <- max( abs( c( sigma - sigma0 )) )
		cat("|     max. parm. change" , round( a1s , 5),"\n")
	
		
        #************************************
        # estimation of error correlation parameters
		a1e <- 0
		epsG <- NULL
		if ( est.corrs ){
			eps.corr0 <- eps.corr
			# identify different eps parameter groups
			epsG <- setdiff( unique( est.eps.corr ) , 0 )
			prbar <- seq( 1 , 10 , len = length(epsG) )
			prbar <- floor( prbar )
			prbar <- c( prbar[1] , diff(prbar) )
			cat(" Estimation of eps:   |")
			for (ee in epsG){
				# bb <- bG[1]
	#            est.bb <- est.b * (est.b == bb )
				est.ee <- 1*(est.eps.corr == ee )
					respml0 <- .ll.rasch.pml.probit( b , a , sigma , 
								Q , eps.corr , itempairs  , IP , eps=10^(-14) )
					ll0 <- respml0$ll
					respml <- .ll.rasch.pml.probit( b , a , sigma , 
								Q,eps.corr + h*est.ee , itempairs  , IP , eps=10^(-14) )
					 ll1 <- respml$ll
					 respml <- .ll.rasch.pml.probit( b , a , sigma , 
							Q , eps.corr - h*est.ee , itempairs  , IP , eps=10^(-14) )
					 ll2 <- respml$ll                           
				d1 <- ( ll1 - ll2  ) / ( 2 * h )    
				# second order derivative
				# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
				d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2               
				if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
				eps.change <- - d1 / d2
				eps.change <- ifelse( abs( eps.change ) > .35 , .35*sign(eps.change) , eps.change )             
				eps.change <- eps.change * est.ee
				eps.corr <- eps.corr + eps.change
				eps.corr[ eps.corr < h & est.eps.corr != 0] <- 2*h
				eps.corr[ eps.corr > 1-h & est.eps.corr != 0 ] <- 1 - 2*h
				cat( paste( rep( "-" , prbar[match(ee,epsG)]), collapse="") )
				flush.console()             
						}
			a1e <- max( abs( eps.corr - eps.corr0 ) )
			cat("|     max. parm. change" , round( a1e , 5),"\n")
						}
		######################################
        # convergence display
        a1 <- stats::aggregate( b , list( est.b) , mean )
        a1aa <- stats::aggregate( a , list( est.a) , mean )		
#        cat("   b parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )
#        cat("   a parameters: " , paste( round( a1aa[,2] , 3 ) , collapse= " " ) , "\n" )
#		if ( D== 1){ 		
#			cat("   sigma parameter:  " , paste( round( sigma, 3 ) , collapse= " " ) , "\n" )
#				}
		if (D>1){
			cat("   sigma parameters:  " , paste( round( sigma[ ! lower.tri(sigma) ] , 3 ) , collapse= " " ) , "\n" )
					}
#		if (est.corrs){
#			a1 <- aggregate( eps.corr , list( est.eps.corr) , mean )
#			cat("   eps parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )		
#					}
        #******************************************************************************
        iter <- iter + 1 
        dev <- -2*ll0
        dev.change <- abs( ( dev - dev0)/ dev0 )
#        par.change <- max( a1a , a1b , a1d , a1k , a1m , a1s)
        par.change <- max( a1b , a1s , a1e , a1a )
        cat( "Pseudolikelihood Objective Function = "  ,   round( dev , 5 ) , "| max. parm. change = " , 
                                        round( par.change , 6 ) ,  " \n"   )  
        if ( ( dev > dev0 ) & ( iter > 4 ) ){ cat("   Objective Function has increased! Convergence Problems?\n") }
        flush.console()        
                }
		#********************************************************
	# information criteria
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat0) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
        ic[[ "np" ]] <- length(bG) + est.sigma + length( epsG )
        # AIC
   #     ic$AIC <- dev + 2*ic$np
        # BIC
    #    ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC
     #   ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np	
	#**********************************************************************************
	# results item parameters			
	item <- data.frame( "item" = colnames(dat0) , 
				"N" = colSums(!is.na(dat0)) , 
				"sumWeights" = colSums( ( !is.na(dat0)) * WM ) , 
				"p" = colMeans( dat0 , na.rm=TRUE ), 
				"b" = b , "est.b"= est.b , 
				"a" = a , "est.a" = est.a )
	if (D==1){ item$sigma <- sigma }							
	item$est.sigma = 1*est.sigma
#	item$link" = link
	if ( ! is.null( itemcluster) ){ item$itemcluster <- itemcluster }
	item$b.logit <- item$b * 1.701
	item$a.logit <- item$a*1.701
	if (D==1){	item$sigma.logit <- item$sigma * 1.701 }
	# add results dependency parameter for item clusters
	#item$itemcluster <- itemcluster
	#item$delta <- 0
    cat("---------------------------------------------------------------------------------------------------------- \n")
	# print item summary
	cat("Item Parameter Summary\n")
	cat( " Estimated" , length(bG) , "Item Parameters\n\n")
	.pr( item , digits=3 )		# print item statistics
    cat("---------------------------------------------------------------------------------------------------------- \n")
 
	#....................................................................
	# print Trait parameter summary
#	cat("Trait Distribution Summary\n")
#	cat("\nCovariance Matrix (Probit Link)\n")
	cat("Trait SD (Probit Link): ")
    cat( round(sigma , 3 ) , "\n")
	cat("Trait SD (Logit Link) : ")
    cat( round( item$sigma[1] * 1.701 , 3 ) , "\n")	
	if (D>1){
		cat("\nCorrelation Matrix\n")
		print( cov2cor(sigma) , digits=3 )
			}
   cat("---------------------------------------------------------------------------------------------------------- \n")
	# print correlation parameter summary
	error.corr0 <- error.corr
	if ( est.corrs){
		cat("Residual Correlation Parameter Summary\n")
		cat( " Estimated" , length(epsG) , "Residual Correlation Parameters\n\n")
		diag(error.corr0) <- 0
		colnames(error.corr) <- rownames(error.corr) <- colnames(dat)
		error.corr[ upper.tri( error.corr) ] <- eps.corr
		error.corr[ lower.tri( error.corr) ] <- 0
		error.corr <- error.corr + t(error.corr)
#		error.corr[ lower.tri( error.corr) ] <- eps.corr
		diag(error.corr) <- 1
	   .pr( error.corr , digits=3 )		# print item statistics
       cat("---------------------------------------------------------------------------------------------------------- \n")
			}
        # computational time
        s2 <- Sys.time()
        if (progress){ 
                cat("---------------------------------------------------------------------------------------------------------- \n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("---------------------------------------------------------------------------------------------------------- \n")
                    }  
	res <- list( "item" = item , "iter" = iter , "deviance" = dev ,
					"b" = b , "sigma" = sigma , 
					"dat" = dat	,  "ic" = ic	, 
					"link" =link , "itempairs" = itempairs ,
					"error.corr" = error.corr0 		,
					"bG" = bG , "aG" = aG , "epsG" = epsG , "est.corrs" = est.corrs ,
					"Q"=Q , "D"=D)	
	if (est.corrs){ 
			res$eps.corr <- eps.corr 
			res$eps.corrM <- error.corr
				}
	class(res) <- "rasch.pml"
	return(res)					
               }
########################################################################## 
                
