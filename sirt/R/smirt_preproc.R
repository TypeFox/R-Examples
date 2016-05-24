
################################################################
	# define estimation functions depending on condensation type
	#-#-#  calculation of probabilities
	.smirt.calcprob <- function( a , b, Q, thetak , c , d , mu.i , irtmodel ){
			if ( irtmodel=="noncomp"){
				res <- calcprob.noncomp( a , b, Q, thetak , c , d )
									}
			if ( irtmodel=="comp"){
				res <- calcprob.comp( a , b, Q, thetak , c , d )
								}
			if ( irtmodel=="partcomp"){
				res <- calcprob.partcomp( a , b, Q, thetak , c , d , mu.i)
								}								
			return(res)
								}
################################################################								
	#-#-#	estimation of b parameters							 
	.smirt.est.b <- function(   b , a , c , d , mu.i , Qmatrix , est.b , theta.k , 
			n.ik , I , K , TP , D ,  numdiff.parm, 
			max.increment ,	msteps ,  mstepconv , irtmodel , increment.factor ){
			if ( irtmodel=="noncomp"){			
			res <-	.smirt.est.b.noncomp(   b , a , c , d , Qmatrix , est.b , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment ,	msteps ,  mstepconv , increment.factor)	
								}
			if ( irtmodel=="comp"){			
			res <-	.smirt.est.b.comp(   b , a , c , d , Qmatrix , est.b , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment ,	msteps ,  mstepconv , increment.factor)	
								}														
			if ( irtmodel=="partcomp"){			
			res <-	.smirt.est.b.partcomp(   b , a , c , d , mu.i , Qmatrix , est.b , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment ,	msteps ,  mstepconv , increment.factor)	
								}									
			return(res)
									}
################################################################
	#-#-#   estimation of a parameters								
	.smirt.est.a <- function(  b , a , c , d , mu.i , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , irtmodel  , increment.factor){		
			if ( irtmodel=="noncomp"){				
			 res <- .smirt.est.a.noncomp(  b , a , c , d , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , increment.factor)		
									}
			if ( irtmodel=="comp"){				
			 res <- .smirt.est.a.comp(  b , a , c , d , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , increment.factor)		
									}	
			if ( irtmodel=="partcomp"){				
			 res <- .smirt.est.a.partcomp(  b , a , c , d , mu.i,  Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , increment.factor)		
									}										
			return(res)
						}
						
################################################################						
	#-#-#   estimation of c parameters										
	.smirt.est.c <- function(  b , a , c , d , mu.i ,  Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , irtmodel  , increment.factor){		
			if ( irtmodel=="noncomp"){				
			 res <- .smirt.est.c.noncomp(  b , a , c , d , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
								}
			if ( irtmodel=="comp"){				
			 res <- .smirt.est.c.comp(  b , a , c , d , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
								}															
			if ( irtmodel=="partcomp"){				
			 res <- .smirt.est.c.partcomp(  b , a , c , d , mu.i , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
								}																	
			return(res)
						}
						
################################################################						
	#-#-#   estimation of d parameters										
	.smirt.est.d <- function(  b , a , c , d , mu.i , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , irtmodel , increment.factor ){		
			if ( irtmodel=="noncomp"){
			 res <- .smirt.est.d.noncomp(  b , a , c , d , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
									}
			if ( irtmodel=="comp"){
			 res <- .smirt.est.d.comp(  b , a , c , d , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
									}
			if ( irtmodel=="partcomp"){
			 res <- .smirt.est.d.partcomp(  b , a , c , d , mu.i , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
									}									
			return(res)
						}
						
						
						
						
########################################################################
.smirt.check.inits <- function( a.init , b.init , irtmodel , Qmatrix){
      #****** check b.init
		if ( ! is.null( b.init) ){
			if (irtmodel=="comp"){
				if ( is.matrix(b.init) ){
					b.init <- b.init[,1]
					cat("*** I only use the first column of the b.init matrix \n***" ,
						 "because a compensatory model is requested!\n" )	
								}
				}  # end comp
			if (irtmodel=="noncomp"){
				if ( is.vector(b.init) ){
					stop("*** A matrix for b.init is required.")
								}
				}  # end noncomp								
			} # end b.init
		#**************
      #****** check b.init
		if ( ! is.null( a.init) ){
				if ( is.vector(a.init) ){
					stop("*** A matrix for a.init is required.")
								}
				if ( is.matrix(a.init) ){
				   if ( ncol(a.init) != ncol(Qmatrix) ){
					stop("*** Check number of dimensions for specifying a.init.")
										}
								}																
				}  # end a.init							
		
	res <- list("b.init"=b.init , "a.init"=a.init)
    return(res)  	
		}
##########################################################



#############################################
# restrict maximum increment
.adj.maxincrement.parameter <- function( oldparm , newparm , 
		max.increment ){
	ISMATR <- is.matrix( oldparm )
	max.increment0 <- max.increment
	if (ISMATR){
		D <- ncol(newparm)	
		for (dd in 1:D){	
	    if ( is.matrix(max.increment) ){
			max.increment <- max.increment0[1,dd] 
						}		
		#	dd <- 1
			increment <- newparm[,dd] - oldparm[,dd]
			increment2 <- ifelse( abs(increment) > max.increment , 
								sign(increment)*max.increment , increment )
			newparm[,dd] <- oldparm[,dd] + increment2
						}
				}
	if (!ISMATR){
			increment <- newparm - oldparm
			increment2 <- ifelse( abs(increment) > max.increment , 
								sign(increment)*max.increment , increment )
			newparm <- oldparm + increment2
						}
    return(newparm)
		}
###########################################################