



#*******************************************************
# Summary for rasch.mml object                         *
##NS S3method(summary,rasch.mml)
summary.rasch.mml <- function( object , ... ){
    # object      ... object from rasch.mml                #
	
	npirt <- object$irtmodel == "npirt"	
	D <- object$D
	
	a5 <- 1*( npirt & ncol( object$item) == 5 )

	cat("------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")

	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")			
		
    cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
	if ( object$Rfcttype == "rasch.mml" ){ cat("Function 'rasch.mml' \n\n") }
	if ( object$Rfcttype == "rasch.mml2" ){ cat("Function 'rasch.mml2' \n\n") }	
	if ( ! object$ramsay.qm ){
		cat("Rasch Type Model with Fixed Discrimination, Guessing and Slipping Parameters \n") 
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		moments <- genlogis.moments( alpha1=object$alpha1 , alpha2=object$alpha2)
		cat("Moments: \n" ); print(round(moments,2)) ; cat("\n")
							}
	    if ( object$ramsay.qm ){  cat("Quotient Model (Ramsay, 1989) \n") 	}
	    if ( object$irtmodel == "npirt" ){  cat("Nonparametric IRT \n") 	}
		if ( object$irtmodel == "missing1"){ cat("Missing Data IRT Model \n")		}
        utils::flush.console()
	if ( sum(object$est.c) > 0){ cat(paste( "Estimated guessing parameter groups \n") )}  
				## estimated guessing parameters
    if ( object$G > 1 ){ cat("\nMultiple Group Estmation with",object$G , "Groups \n") 
			print(object$groupindex) ; cat("\n")
			}
	cat("------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , 
				round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , 
					round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , 
				round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

	#------------------------------------
	if ( object$D == 1){	
	if ( is.null( object$trait.weights) ){
		cat( "Trait Distribution (" , length(object$trait.distr[,1]) , " Knots )\n" , 
				  "Mean=" , round( object$mean.trait,3) , "\n SD=" , round( object$sd.trait , 3) ,
				  "\n Skewness=" , round( object$skewness.trait , 3) 
				  ) 
						}
	if ( ! is.null( object$trait.weights) ){
		M1 <- stats::weighted.mean( object$theta.k ,object$trait.weights )
		S1 <- sqrt( stats::weighted.mean( object$theta.k^2 ,object$trait.weights ) - M1^2	)
		cat( "Fixed Trait Distribution (" , length(object$trait.distr[,1]) , " Knots )\n" , 
				  "Mean=" , round( M1 ,3 ) , 
					" SD=" , round( S1 , 3) 
								)
						}						
						}
	if ( object$ramsay.qm ){ cat( "      Note: log theta distribution is parametrized!") }
	cat("\n")
	if ( D > 1){
		cat("Mean Vector\n") ; print( round( object$mu , 3 ) )
		cat("\nCovariance Matrix\n") ; print( round( object$Sigma.cov , 3 ) )	
		cat("\n")
		covmat <- object$Sigma.cov 
		covmat2 <- stats::cov2cor( covmat )
		diag(covmat2) <- sqrt( diag(covmat) )
		cat("\nStandard Deviations / Correlation Matrix\n") ; print( round( covmat2 , 3 ) )	
		cat("\n")
				}
	
	if ( object$irtmodel != "npirt" ){	
		cat( "Item Difficulty Distribution (" , nrow(object$item) , " Items )\n" , 
				  "Mean=" , round( mean(object$item$b) ,3) , " SD=" , 
							round( stats::sd(object$item$b) , 3) , "\n") 
							}
    cat( "Distribution of Items Administered (" , nrow(object$item) , " Items )\n" , 
              "Mean=" , round( mean(rowSums( 1 - is.na(object$dat) )) ,3) , " SD=" , 
                        round( stats::sd(rowSums( 1 - is.na(object$dat) )) ,3) , "\n\n") 
	cat( "EAP Reliability: ") 
	cat(round( object$reliability$eap.reliability,3 ) )
	cat( "\n")
	if ( a5 == 0 ){
	cat("------------------------------------------------------------\n")
		cat("Item Parameter \n")
		if ( ! object$ramsay.qm ){ obji <- object$item } else { obji <- object$item2 }
		rvars <- seq( 2 , ncol(obji ) )
		ind <- which( colMeans( is.na( obji )) == 1 )
		roundvars <- setdiff( rvars , ind )
		if (npirt ){ 
		  roundvars <- c("p","est","se" )
				}		
		for (vv in roundvars ){ 
			obji[,vv] <- round( obji[,vv] , 3 ) 
					}
		rownames(obji) <- NULL
		print( obji )                
			}
                }
#*******************************************************



