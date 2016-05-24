#*******************************************************
# Summary for rasch.mirtlc object                         *
summary.rasch.mirtlc <- function( object,... ){
    # object      ... object from rasch.mml                #
    cat("---------------------------------------------------------------------------------------------------------- \n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Multidimensional Item Response Latent Class Model \n\n")
	modeltype <- object$modeltype
	if (object$modeltype =="LC" ){
		cat("Latent Class Model with" , object$Nclasses , "Classes" , 
			" - " , object$G , "Group(s)\n")  
					}
	if (object$modeltype %in% c("MLC1","MLC2") ){
		if (modeltype=="MLC1"){
				cat("Multidimensional Latent Class Rasch Model with\n     " )
							}
		if (modeltype=="MLC2"){
				cat("Multidimensional Latent Class 2PL Model with\n     " )
							}
			cat( object$Nclasses , "Classes" , 
				", " , object$G , "Group(s)", "," ,
				object$D , "Dimension(s)\n")  
					}
	if (object$distribution.trait=="normal" ){
		cat("      Normal distribution assumption\n" ) 
					}	
	if (object$distribution.trait=="smooth2" ){
		cat("      Log-linear Smoothing (2 Moments)\n" ) 
					}
	if (object$distribution.trait=="smooth3" ){
		cat("      Log-linear Smoothing (3 Moments)\n" ) 
					}
	if (object$distribution.trait=="smooth4" ){
		cat("      Log-linear Smoothing (4 Moments)\n" ) 
					}
    cat("---------------------------------------------------------------------------------------------------------- \n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    

    cat( "Number of estimated item parameters = " , object$ic$itempars , "\n" )    
    cat( "Number of estimated distribution parameters = " , object$ic$traitpars , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
	
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

		cat("---------------------------------------------------------------------------------------------------------- \n")
		cat("Trait Distribution\n")
		obji <- object$pi.k
		if ( object$modeltype %in% c("MLC1","MLC2") ){ obji <- object$trait }
		obji <- round( obji , 3)
		print( obji )                		
		if ( object$modeltype %in% c("MLC1","MLC2") ){ 
			cat( "\nM Trait: " , round( object$mean.trait , 3 ) , "\n")
			cat( "SD Trait: " , round( object$sd.trait , 3 ) , "\n")
			cat( "Skewness Trait: " , round( object$skewness.trait , 3 ) , "\n")
			if ( object$D > 1){
					cat( "Correlations Trait: \n" )
					for (gg in 1:object$G){
						cat("Group" , gg , "\n")
						print( round( object$cor.trait[,,gg] , 3 ) )
								}
									}
									}
		cat("---------------------------------------------------------------------------------------------------------- \n")
		cat("Item Parameters \n")
		cat("Item Probabilities\n")
		obji <- t(object$pjk)
		obji <- round( obji , 3)
		print( obji )           
		if ( object$modeltype %in% c("MLC1","MLC2") ){
			cat("\nItem Parameter \n")
			obji <- object$item
			obji <- round( obji , 3)
			print( obji )           
					}			
			}
#*******************************************************
