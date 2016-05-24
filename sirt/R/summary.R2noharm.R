
#..........................................................
# R2Noharm summary function

#----------------------------------------------------------
# summary function for R2noharm
summary.R2noharm <- function( object , logfile=NULL , ...){
	#**********************
	# global summary
	if ( ! is.null(logfile) ){ sink( paste0(logfile,".Rout") , split=TRUE ) }
	cat( paste( "-------------------------------------------------") , "\n" )	
	cat( paste( "NOHARM4 " , object$systime , sep="") , "\n" )	
	d <- utils::packageDescription("sirt")
	cat(paste(d$Package," " , d$Version," (",d$Date,")",sep=""))
	cat( paste( "\n-------------------------------------------------") , "\n\n" )	
	#**********************
	# model type
	if (object$modtype == 2){ cat( "Multidimensional Exploratory Factor Analysis\n\n") }
	if (object$modtype == 3){ cat( "Multidimensional Confirmatory Factor Analysis\n\n") }
	if (object$modtype == 4){ cat( "Unidimensional Confirmatory Factor Analysis\n\n") }	
	#*********************
	# descriptives
	cat( paste( "Number of Observations: " , object$Nobs , sep="") , "\n" )	
	cat( paste( "Number of Items       : " , object$Nitems , sep="") , "\n" )		
	cat( paste( "Number of Dimensions  : " , object$dimensions , sep="") , "\n" )		
	cat( paste( "Tanaka Index          : " , round(object$tanaka,object$display.fit) , sep="") , "\n" )
	cat( paste( "RMSR                  : " , round(object$rmsr,object$display.fit) , sep="") , "\n\n" )
	
	# chi square statistic
	if (object$modtype %in% 2:4){
		cat("Chi Square Statistic of Gessaroli & De Champlain (1996)\n\n")
		cat( paste( "Number of Estimated Parameters : " , object$Nestpars , sep="") , "\n" )			
		cat( paste( "Chi2                           : " , round(object$chisquare,3) , sep="") , "\n" )			
		cat( paste( "Degrees of Freedom (df)        : " , object$df , sep="") , "\n" )		
		cat( paste( "p(Chi2,df)                     : " , round(object$p.chisquare,3) , sep="") , "\n" )			
		cat( paste( "Chi2 / df                      : " , round(object$chisquare / object$df,3) , sep="") , "\n" )			
		cat( paste( "RMSEA                          : " , round(object$rmsea,3) , sep="") , "\n\n" )				
						}
	if ( object$modtype %in% 2:4){ 
		#***
		# factor correlation
		cat( "Factor Covariance Matrix\n")				
		print( object$factor.cor )
		if ( object$modtype %in% 3){ 
			cat( "\nFactor Correlation Matrix\n")				
			print( round( stats::cov2cor(object$factor.cor) , 3 ))
									}
								}
	if ( object$modtype %in% 3:4){
		# item parameter
		l1 <- object$loadings.theta
		cat("\nItem Parameters - Latent Trait Model Parametrization\n",
				"Loadings, Constants, Guessing Parameters and Decriptives\n\n")
		m1 <- l1 %*% as.matrix(object$factor.cor) %*% t( l1 )
		v1 <- 1 + diag( m1 )
		l1 <- data.frame( l1 , "final.constants" = object$final.constants , "guess" = object$guesses ,
					"item.variance" = round( v1 , 3) , "N" = diag(object$N.itempair) , "p" = diag(object$pm) )
		print( l1 )
		res <- list( "itempars" = l1 )
		# factor analysis parametrization
		l1 <- object$loadings
		cat("\nItem Parameters - Common Factor Parametrization\n",
			"Loadings, Thresholds, Uniquenesses and Guessing Parameters\n\n")
		l1 <- data.frame( l1 , "threshold" = object$thresholds , "uniqueness" = object$uniquenesses , 
						"guess" = object$guesses )		
		print( l1 )
		res[["itempars.faparm" ]] <- l1 
		#***
		# Lord's parametrization for unidimensional models
		if (object$modtype == 4){
			cat("\nItem Parameters - Lord's Parametrization\n",
				"Difficulty and Discrimination\n\n")		
			l1 <- data.frame( "difficulty" = object$difficulties , 
							   "discrimination" = object$discriminations )
			print(l1)
		res[["itempars.lord" ]] <- l1 
					}
			}
	if ( object$modtype %in% 2){
		# item parameter
		l1 <- object$loadings.theta
		cat("\nItem Parameters - Promax Rotated Parameters\n",
				"Loadings, Constants, Guessing Parameters and Descriptives\n\n")
		l1 <- object$promax.theta
		l1 <- data.frame( l1 , "final.constants" = object$final.constants , 
					"guess" = object$guesses , "N" = diag(object$N.itempair) , "p" = diag(object$pm) )
		print( l1 )
		res <- list( "itempars.promax" = l1 )
				}
	if ( ! is.null(logfile) ){ sink() }				
		invisible(res)
	}
#--------------------------------------------------------------------------------------
