

##################################################################
# Summary of the GDINA model
summary.gdina <- function( object , digits = 4 , file = NULL , ... ){
	#-------------------------------------------------------
	# INPUT:
	# object	... result from GDINA analysis
	# rdigits 	... number of digits for rounding parameter estimates
	#-------------------------------------------------------
	rdigits <- digits
	
 	osink( file = file , suffix = paste0( file, "__SUMMARY.Rout") )


	# Parameter summary
    cat("---------------------------------------------------------------------------------------------------------- \n")
	d1 <- utils::packageDescription("CDM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( "Date of Analysis:" , paste( object$time$s2 ) , "\n" )
	cat("Computation Time:" , print(object$time$s2 - object$time$s1), "\n\n")	
    if (object$HOGDINA==-1){ 
		cat("Generalized DINA Model \n") } else {
		cat("Higher Order Generalized DINA Model \n") }
    if ( object$G > 1 ){ 
			cat("  Multiple Group Estmation with",object$G , "Groups \n") 
			# group statistics
			cat("\nGroup statistics\n")
			print( object$group.stat )				
			cat("\n")
						}
						
	cat( "\nNumber of iterations =" , object$iter  )
	if ( ! object$converged ){ cat("\nMaximum number of iterations was reached.\n") }		
	cat( "\nIteration with minimal deviance =" , object$iter , "\n" )	
	
	
    cat( "\nDeviance =" , round( object$deviance , 2 ) ) 
	cat( "  | Loglikelihood =" , round( - object$deviance / 2 , 2 ) ,	"\n" )
    cat( "Number of persons =" , object$N , "\n" )    
	    cat( "Number of items =" , ncol(object$data) , "\n" )    
    cat( "Number of estimated parameters =" , object$Npars , "\n" )    
    cat( "Number of estimated item parameters =" , object$Nipar , "\n" )    	
    cat( "Number of estimated skill class parameters =" , object$Nskillpar )    		
	cat( " (" , object$Nskillclasses , "latent skill classes)\n")
    cat( "\nAIC = " , round( object$AIC , 2 ) , " ; penalty =" , round( object$AIC - object$deviance ,2 ) , "\n" )    
    cat( "BIC = " , round( object$BIC , 2 ) , " ; penalty =" , round( object$BIC - object$deviance ,2 ) , "\n" )  
    cat( "CAIC = " , round( object$CAIC , 2 ) ," ; penalty =" , round( object$CAIC - object$deviance ,2 ) , "\n\n" )  
#	if (object$reduced.skillspace ){ 
#	    cat("Goodness of fit for reduced skillspace\n")
#		cat( "Delta index =" , round( object$delta.index , 3 ) , "\n\n")
#						}

#    cat("Model fit\n")
#	g1 <- gdina.fit( object , print.output = TRUE )				
	###########################################################
	ds <- object$coef
	selvars <- intersect( c("est", "se" ) , colnames(ds) )
	ind <- which( colnames(ds) %in% selvars ) 
#	if (G>0){ ind <- which( colnames(ds) %in% c("est" ) ) }
    for (ii in ind){
		ds[,ii] <- round( ds[,ii] , rdigits )
					}
	cat("----------------------------------------------------------------------------\n")
	cat("\nItem Parameter Estimates \n\n")
	r1 <- options()
	options(scipen=999)
	print(ds)
	options(scipen=r1$scipen)
    	if ( ! is.null( object$delta.designmatrix ) ){ 
			cat("\nNote: Standard errors are not (yet) correctly implemented!\n")
											}	

	cat("\nRMSEA Item Fit\n")											
	print( round( object$itemfit.rmsea,3) )
											
    cat("\nMean of RMSEA item fit:" , 
     round( object$mean.rmsea ,3 ) , "\n")											
	
	# RRUM model	
	if (object$rrum.model){
		cat("\n****\nRRUM Parametrization\n")
		print( round( object$rrum.params,3) , na  ="")
		cat("\n")
					}

	cat("----------------------------------------------------------------------------\n")
	cat("Model Implied Conditional Item Probabilities \n\n")
	obji <- object$probitem	
	obji[,"prob"] <- round( obji$prob, 4 )	
	print(obji)					
	cat("----------------------------------------------------------------------------\n")
	cat("\nSkill Probabilities \n\n")
	print(round(object$skill.patt ,rdigits) )
	
	#**** output tetrachoric or polychoric correlations
		cat("----------------------------------------------------------------------------\n")

		cat("\nPolychoric Correlations \n")
		G <- object$G
		for (gg in 1:G){
			cat( paste0( "\nGroup ", gg , "\n") )
			obji <- object$polychor[[gg]]	
			print( round( obji , 3 ))
						}

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

		if (object$HOGDINA>=0){
			cat("\n***************************\n")	
			cat("Higher Order GDINA Model ")
			cat("\n  Attribute Response Function Parameters \n\n")
			print( round( object$attr.rf,3) )
			}

   csink( file = file )	
			
		}
##########################################################################

