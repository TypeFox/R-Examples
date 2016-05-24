
############################################
# probability in compensatory model
## extern "C" {
## SEXP SMIRT_CALCPROB_comp( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
calcprob.comp <- function (a,b,Q,thetak,cc,dd){ 
	.Call("SMIRT_CALCPROB_COMP", a,b,Q,thetak,cc,dd, PACKAGE = "sirt")
					}


###########################################
# estimation of b
.smirt.est.b.comp <- function(   b , a , c , d , Qmatrix , est.b , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=1,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.b		# zeros are allowed!
	cat("  M steps b parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- rep(1,I)
	Q2[ est.b == 0 ] <- 0
	Q <- Qmatrix	
	b00 <- b
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		b0 <- b

			probres <- calcprob.comp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a , b + h*Q2, Q , thetak=theta.k , c , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a , b - h*Q2, Q , thetak=theta.k , c , d )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))	)
			b <- b + (res$increment)[ind]
			se.b <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( b - b0 ) )
		it <- it+1
		cat("-") 
			}
	cat(" " , it , "Step(s) \n")	
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
.smirt.est.a.comp <- function(   b , a , c , d , Qmatrix , est.a , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.a.increment=max.a.increment,
		msteps ,  mstepconv  , increment.factor){
    h <- numdiff.parm
	diffindex <- est.a		# zeros are allowed!
	Q <- Qmatrix
	cat("  M steps a parameter   |")
	it <- 0 ;	conv1 <- 1000	
	Q2 <- Q1 <- 0*Qmatrix
	se.a <- a
	a00 <- a
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a0 <- a
		
		for (dd in 1:D){
#				dd <- 2
			Q2 <- Q1
			Q2[,dd] <- 1 * ( Qmatrix[,dd] != 0 )
			
			probres <- calcprob.comp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a + h*Q2 , b , Q , thetak=theta.k , c , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a - h*Q2 , b , Q , thetak=theta.k , c , d )
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
		cat("-") 
			}
	if ( increment.factor > 1){
		a <- .adj.maxincrement.parameter( oldparm=a00 , newparm=a , 
					max.increment=max.a.increment )		
						}				
	cat(" " , it , "Step(s) \n")	
    res <- list("a" = a , "se.a" = se.a , 
			"ll" = sum(res$ll0) )
    return(res)
			}				


###########################################
# estimation of c
.smirt.est.c.comp <- function(   b , a , c , d , Qmatrix , est.c , theta.k , 
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
			probres <- calcprob.comp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a  , b , Q , thetak=theta.k , c+h*Q1 , d )
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a  , b , Q , thetak=theta.k , c-h*Q1 , d )
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
		cat("-") 
			}
	cat(" " , it , "Step(s) \n")	
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
.smirt.est.d.comp <- function(   b , a , c , d , Qmatrix , est.d , theta.k , 
        n.ik , I , K , TP , D , numdiff.parm=.001 , max.increment=max.increment,
		msteps ,  mstepconv  , increment.factor){
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
			probres <- calcprob.comp( a , b, Q , thetak=theta.k , c , d )
			pjk <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a  , b , Q , thetak=theta.k , c , d +h*Q1)
			pjk1 <- problong2probarray( probres , I , TP )		

			probres <- calcprob.comp( a  , b , Q , thetak=theta.k , c , d-h*Q1 )
			pjk2 <- problong2probarray( probres , I , TP )		
			
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.increment , numdiff.parm )
			ind <- match( diffindex , sort(unique( diffindex ))	)
			incr <- res$increment
			incr[ is.na(incr ) ] <- 0			
			d <- d + incr[ind]
			d[ d>1 ] <- .999
#			max.a.increment[,dd] <- abs( (res$increment)[ind] )
			se.d <- (sqrt(  1 / abs(res$d2) ))[ind]	
		conv1 <- max( abs( d - d0 ) )
		it <- it+1
		cat("-") 
			}
	cat(" " , it , "Step(s) \n")	
	if ( increment.factor > 1){
		d <- .adj.maxincrement.parameter( oldparm=d00 , newparm=d , 
					max.increment=max.increment )		
						}		
    res <- list("d" = d , "se.d" = se.d , 
			"ll" = sum(res$ll0) )
    return(res)
			}				
			
