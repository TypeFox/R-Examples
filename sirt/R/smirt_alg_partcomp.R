
############################################
# probability in noncompensatory model
## extern "C" {
## SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
calcprob.partcomp <- function (a,b,Q,thetak,cc,dd,mu.i){ 
	#     Rcpp::NumericMatrix A(a);  
	#     Rcpp::NumericMatrix B(b);  
	#     Rcpp::NumericMatrix QQ(Q);  
	#     Rcpp::NumericMatrix THETA(thetak);  
	#     Rcpp::NumericVector CC(cc);  
	#     Rcpp::NumericVector DD(dd);  
	#     Rcpp::NumericVector MUI(mui) ;
	#**********
	# SEXP SMIRT_CALCPROB_PARTCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ,
	#	SEXP mui ){
	.Call("SMIRT_CALCPROB_PARTCOMP", 
			a,b,Q,thetak,cc,dd, mu.i , 
			PACKAGE = "sirt")
					}


		
###########################################
# estimation of b
.smirt.est.b.partcomp <- function(   b , a , c , d , mu.i , Qmatrix , est.b , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=1,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.b		# zeros are allowed!
	cat("  M steps b parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- Q1 <- 0*Qmatrix
	Q <- Qmatrix	
	se.b <- b
	b00 <- b
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		b0 <- b
		for (dd in 1:D){
#				dd <- 2
			Q2 <- Q1
			Q2[,dd] <- 1 * ( Qmatrix[,dd] != 0 )
			
			probres <- calcprob.partcomp( a , b, Q , thetak=theta.k , c , d , mu.i)
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a , b + h*Q2, Q , thetak=theta.k , c , d , mu.i)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a , b - h*Q2, Q , thetak=theta.k , c , d , mu.i)
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex[,dd] , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex[,dd] , sort(unique( diffindex[,dd] )) )
			b[,dd] <- b[,dd] + (res$increment)[ind]
			se.b[,dd] <- (sqrt(  1 / abs(res$d2) ))[ind]	
						}   # end dd
		conv1 <- max( abs( b - b0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		b <- .adj.maxincrement.parameter( oldparm=b00 , newparm=b , 
					max.increment=max.increment )		
						}
    res <- list("b" = b , "se.b" = se.b , 
			"ll" = sum(res$ll0) )
    return(res)
			}			
		
###########################################
# estimation of a
.smirt.est.a.partcomp <- function(   b , a , c , d ,  mu.i ,Qmatrix , est.a , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.a.increment ,
		msteps ,  mstepconv , increment.factor){
    h <- numdiff.parm
	diffindex <- est.a		# zeros are allowed!
	cat("  M steps a parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- Q1 <- 0*Qmatrix
	Q <- Qmatrix	
	se.a <- a
	a00 <- a
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a0 <- a		
		for (dd in 1:D){
#				dd <- 2
			Q2 <- Q1
			Q2[,dd] <- 1 * ( Qmatrix[,dd] != 0 )
			
			probres <- calcprob.partcomp( a , b, Q , thetak=theta.k , c , d , mu.i )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a + h*Q2 , b , Q , thetak=theta.k , c , d , mu.i)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a - h*Q2 , b , Q , thetak=theta.k , c , d , mu.i)
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex[,dd] , 
					max.increment=max.a.increment[,dd] , numdiff.parm )
			ind <- match( diffindex[,dd] , sort(unique( diffindex[,dd] )) )
			a[,dd] <- a[,dd] + (res$increment)[ind]
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.a[,dd] <- (sqrt(  1 / abs(res$d2) ))[ind]			
						}   # end dd
		conv1 <- max( abs( a - a0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	#****
	# post-processing of a parameters	
	if ( increment.factor > 1){
		a <- .adj.maxincrement.parameter( oldparm=a00 , newparm=a , 
					max.increment=max.a.increment )		
						}	
    res <- list("a" = a , "se.a" = se.a , "ll" = sum(res$ll0) )
    return(res)
			}				
		
			
			
###########################################
# estimation of c
.smirt.est.c.partcomp <- function(   b , a , c , d , mu.i , Qmatrix , est.c , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.c		# zeros are allowed!
	Q1 <- rep(1,I)
	Q1[ est.c == 0 ] <- 0
	Q <- Qmatrix	
	c00 <- c
	cat("  M steps c parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		c0 <- c
			probres <- calcprob.partcomp( a , b, Q , thetak=theta.k , c , d , mu.i)
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c+h*Q1 , d , mu.i)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c-h*Q1 , d , mu.i)
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))	)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0
			c <- c + incr[ind]
			c[ c< 0 ] <- .001
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.c <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( c - c0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		c <- .adj.maxincrement.parameter( oldparm=c00 , newparm=c , 
					max.increment=max.increment )		
						}	
    res <- list("c" = c , "se.c" = se.c , 
			"ll" = sum(res$ll0) )
    return(res)
			}				
			

###########################################
# estimation of c
.smirt.est.d.partcomp <- function(   b , a , c , d , mu.i , Qmatrix , est.d , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv , increment.factor){
    h <- numdiff.parm
	diffindex <- est.d		# zeros are allowed!
	Q1 <- rep(1,I)
	Q1[ est.d == 0 ] <- 0
	Q <- Qmatrix	
	d00 <- d
	cat("  M steps d parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		d0 <- d
			probres <- calcprob.partcomp( a , b, Q , thetak=theta.k , c , d , mu.i)
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c , d +h*Q1 , mu.i)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c , d-h*Q1 , mu.i)
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))		)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0			
			d <- d + incr[ind]
			d[ d>1 ] <- .999
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.d <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( d - d0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		d <- .adj.maxincrement.parameter( oldparm=d00 , newparm=d , 
					max.increment=max.increment )		
						}		
    res <- list("d" = d , "se.d" = se.d , "ll" = sum(res$ll0) )
    return(res)
			}				
			
			
			
	
###########################################
# estimation of mu.i
.smirt.est.mu.i.partcomp <- function(   b , a , c , d , mu.i , Qmatrix , est.mu.i , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.mu.i		# zeros are allowed!
	Q1 <- rep(1,I)
	Q1[ est.mu.i == 0 ] <- 0
	Q <- Qmatrix	
	c00 <- mu.i
	cat("  M steps mu.i parameter   |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		c0 <- mu.i
			probres <- calcprob.partcomp( a , b, Q , thetak=theta.k , c , d , mu.i)
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c , d , mu.i+h*Q1)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.partcomp( a  , b , Q , thetak=theta.k , c , d , mu.i-h*Q1)
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))	)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0
			mu.i <- mu.i + incr[ind]
			mu.i[  mu.i < 0 ] <- .001
			mu.i[  mu.i > 1 ] <- 1 - .001			
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.c <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( mu.i - c0 ) )
		it <- it+1
		cat("-") # ; flush.console()
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
	if ( increment.factor > 1){
		mu.i <- .adj.maxincrement.parameter( oldparm=c00 , newparm= mu.i  , 
					max.increment=max.increment )		
						}	
    res <- list("mu.i" = mu.i , "se.c" = se.c , 
			"ll" = sum(res$ll0) )
    return(res)
			}				
			