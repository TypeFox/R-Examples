

##################################################################
# Summary of the GDINA model
summary.mcdina <- function( object , digits = 4 , file = NULL , ... ){
	#-------------------------------------------------------
	# INPUT:
	# object	... result from GDINA analysis
	# rdigits 	... number of digits for rounding parameter estimates
	#-------------------------------------------------------
	rdigits <- digits
	
 	osink( file = file , suffix = paste0( "__SUMMARY.Rout") )


	# Parameter summary
	cat("----------------------------------------------------------------------------\n")
	d1 <- utils::packageDescription("CDM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( "Date of Analysis:" , paste( object$time$s2 ) , "\n" )
	cat("Computation Time:" , print(object$time$s2 - object$time$s1), "\n\n")	
 	cat("Multiple Choice DINA Model (MC-DINA)\n") 
	cat( "\nNumber of iterations =" , object$iter )
	if ( ! object$converged ){ cat("\nMaximum number of iterations was reached.") }
	
	cat( "\n\nDeviance =" , round( object$ic$dev , 2 ) ) 
	cat( "  | Loglikelihood =" , round( - object$ic$dev / 2 , 2 ) ,	"\n" )
    cat( "Number of persons =" , object$ic$n , "\n" )    
    cat( "Number of groups =" , object$G , "\n" )    	
	    cat( "Number of items =" , ncol(object$dat) , "\n" )    
    cat( "Number of estimated parameters =" , object$ic$np , "\n" )    
    cat( "Number of estimated item parameters =" , object$ic$itempars , "\n" )    	
    cat( "Number of estimated skill class parameters =" , object$ic$traitpars )    		
	cat( " (" , object$ic$Nskillclasses , "latent skill classes)\n")
    cat( "\nAIC = " , round( object$ic$AIC , 2 ) , " ; penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
    cat( "BIC = " , round( object$ic$BIC , 2 ) , " ; penalty =" , 
				round( object$ic$BIC - object$ic$deviance ,2 ) , "\n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," ; penalty =" , 
				round( object$ic$CAIC - object$ic$deviance ,2 ) , "\n\n" )  
#	if (object$reduced.skillspace ){ 
#	    cat("Goodness of fit for reduced skillspace\n")
#		cat( "Delta index =" , round( object$delta.index , 3 ) , "\n\n")
#						}

#    cat("Model fit\n")
#	g1 <- gdina.fit( object , print.output = TRUE )				
	###########################################################
	ds <- object$item
	cds <- colnames(ds)
	inds <- grep( "se" , cds )
	ds <- ds[ , - inds ]
	ind <- grep( "Cat" , colnames(ds) )
	ds[,ind] <- round( ds[,ind] , rdigits )
	cat("----------------------------------------------------------------------------\n")
	cat("\nItem Parameter Estimates \n\n")
	r1 <- options()
	options(scipen=999)
	print(ds)
	options(scipen=r1$scipen)	
	cat("----------------------------------------------------------------------------\n")
	cat("\nSkill Probabilities \n\n")
	print(round(object$skill.patt ,rdigits) )
#    cat("not yet implemented!\n")
#	if ( ( object$G == 1 ) & (ncol(object$q.matrix ) > 1 ) & 
#			max(object$NAttr ==1 ) ){ 
#		cat("----------------------------------------------------------------------------\n")
#		QM <- max(object$q.matrix)
#		if (QM == 1){ 
#			cat("\nTetrachoric Correlations \n\n")
#			gt1 <- skill.cor( object )
#				} else {
#			cat("\nPolychoric Correlations \n\n")
#			gt1 <- skill.polychor( object )				
#				}
#		print(round(gt1$cor.skills,3))
#			}
	cat("\n----------------------------------------------------------------------------\n")	
	cat("\nSkill Pattern Probabilities \n\n")
	if ( object$G == 1 ){
		xt <- round( object$attribute.patt[,1] , rdigits )
		names(xt) <- rownames( object$attribute.patt )
			} else {
	xt <- round( object$attribute.patt , rdigits )
	rownames(xt) <- rownames( object$attribute.patt )
					}
	print(xt)

   csink( file = file )		
		
		}
##########################################################################


#***************************************************************
# RRUM parametrization
.rrum.param <- function( delta.summary , q.matrix ){
	#---
	#  RRUM parametrization
	#  log( P(X=1) ) = b0 + b1*alpha1 + b2 * alpha2 
	#  RRUM:
	#  P(X=1) = pi * r1^( 1- alpha1) * r2^(1-alpha2)
	#  => log( P(X=1) ) = log[ pi * r1 * r2 * r1^(-alpha1) * r2^(-alpha2) ]
	#                   = log( pi ) + log(r1) + log(r2) + -log(r1)*alpha1 + -log(r2) * alpha2
	#  => b1 = -log(r1) and r1 = exp( -b1 )
	#  => log(pi) = b0 + b1 + b2 and pi = exp( b0 + b1 + b2 )
	I <- nrow(q.matrix)
	K <- ncol(q.matrix)
	rrum.params <- matrix( NA , I , K+1 )
	rownames(rrum.params) <- delta.summary[ delta.summary$partype == 0 , "item" ]
	colnames(rrum.params) <- c( "pi" , paste( "r_", colnames(q.matrix) , sep="") )
	for (ii in 1:I){
		# ii <- 2
		d.ii <- delta.summary[ delta.summary$itemno == ii , ]
		rrum.params[ii,"pi"] <- exp( sum( d.ii$est ) )
		rrum.params[ ii , which( q.matrix[ii,]==1) +1 ] <- exp( - d.ii$est[-1] )
				}
	return( rrum.params )
        }