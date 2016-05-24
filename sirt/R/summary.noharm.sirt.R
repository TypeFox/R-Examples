
#..........................................................
# R2Noharm summary function

#----------------------------------------------------------
# summary function for R2noharm
summary.noharm.sirt <- function( object , logfile=NULL , ...){
	#**********************
	# global summary
	if ( ! is.null(logfile) ){ sink( paste0(logfile,".Rout") , split=TRUE ) }
	cat( paste( "-------------------------------------------------") , "\n" )	
#	cat( paste( "NOHARM | '", class(object),"' | " , object$systime , sep="") , "\n" )	
	cat( paste( "Function '", class(object),"' | " , object$systime , sep="") , "\n" )	
	d <- utils::packageDescription("sirt")
	cat(paste(d$Package," " , d$Version," (",d$Date,")",sep=""))
	cat( paste( "\n-------------------------------------------------") , "\n\n" )	
	#**********************
	# model type
	if (object$modtype == 2){ cat( "Multidimensional Exploratory Factor Analysis\n") }
	if (object$modtype == 3){ cat( "Multidimensional Confirmatory Factor Analysis\n") }

	if (object$modesttype==1){
		cat( "NOHARM approximation\n\n") 
					}
	if (object$modesttype==2){
		cat( "Estimation based on tetrachoric correlations\n\n") 
					}

					
	#*********************
	# descriptives
	cat( paste( "Number of Observations: " , object$Nobs , sep="") , "\n" )	
	cat( paste( "Number of Items       : " , object$Nitems , sep="") , "\n" )		
	cat( paste( "Number of Dimensions  : " , object$dimensions , sep="") , "\n" )		
	cat( paste( "Number of Iterations  : " , object$iter , sep="") , "\n" )	
	if (object$modesttype==1){
		cat( paste( "Tanaka Index          : " , round(object$tanaka,object$display.fit) , sep="") , "\n" )
		cat( paste( "RMSR                  : " , round(object$rmsr,object$display.fit) , sep="") , "\n\n" )
					} else { cat("\n") }

		cat( paste( "Number of Used Item Pairs      : " , object$sumwgtm , sep="") , "\n" )		
		cat( paste( "Number of Estimated Parameters : " , object$Nestpars$total , sep="") , "\n" )			
		cat( paste( "       # Thresholds            : " , object$Nestpars$thresh , sep="") , "\n" )		
		cat( paste( "       # Loadings              : " , object$Nestpars$F , sep="") , "\n" )		
		cat( paste( "       # Variances/Covariances : " , object$Nestpars$P , sep="") , "\n" )				
		cat( paste( "       # Residual Correlations : " , object$Nestpars$Psi , sep="") , "\n\n" )
		

		
	# chi square statistic
	if (object$modtype %in% 2:4){
	   if (object$modesttype==1){
		cat("Chi Square Statistic of Gessaroli & De Champlain (1996)\n\n")				
		cat( paste( "Chi2                           : " , round(object$chisquare,3) , sep="") , "\n" )			
		cat( paste( "Degrees of Freedom (df)        : " , object$df , sep="") , "\n" )		
		cat( paste( "p(Chi2,df)                     : " , round(object$p.chisquare,3) , sep="") , "\n" )			
		cat( paste( "Chi2 / df                      : " , round(object$chisquare / object$df,3) , sep="") , "\n" )			
		cat( paste( "RMSEA                          : " , round(object$rmsea,3) , sep="") , "\n\n" )
				}
		cat( paste( "Green-Yang Reliability Omega Total : " , round(object$omega.rel,3) , 
					sep="") , "\n\n" )		
						}
	if ( object$modtype %in% 2:4){ 
		#***
		# factor correlation
		cat( "Factor Covariance Matrix\n")				
		print( round( object$factor.cor ,3 ) )
		if ( object$modtype %in% 3){ 
			cat( "\nFactor Correlation Matrix\n")				
			print( round( cov2cor(object$factor.cor) , 3 ))
									}
								}
	if ( object$modtype %in% 3:4){
		# item parameter
		l1 <- object$loadings.theta
		cat("\nItem Parameters - Latent Trait Model (THETA) Parametrization\n",
				"Loadings, Constants, Asymptotes and Decriptives\n\n")
		m1 <- l1 %*% as.matrix(object$factor.cor) %*% t( l1 )
		v1 <- 1 + diag( m1 )
		l1 <- data.frame( l1 , "final.constant" = object$final.constants , 
					"lower" = object$lower , "upper" = object$upper ,
					"item.variance" = round( v1 , 3) , "N" = diag(object$N.itempair) , 
					"p" = diag(object$pm) )
		l1 <- round(l1,3)
		print( l1 )
		res <- list( "itempars" = l1 )
		# residual correlation
		if (object$estpars$estPsi==1){
			cat("\nResidual Correlation Matrix\n")		
			l1 <- round( object$residcorr , 3 )
			print(l1)
			}				
		# factor analysis parametrization
		l1 <- object$loadings
		cat("\nItem Parameters - Common Factor (DELTA) Parametrization\n",
			"Loadings, Thresholds, Uniquenesses and Asymptotes\n\n")
		l1 <- data.frame( l1 , "threshold" = object$thresholds , 
						"lower" = object$lower , "upper"=object$upper ,
						"uniqueness" = object$uniquenesses )
		l1 <- round(l1,3)						
		print( l1 )
		res[["itempars.faparm" ]] <- l1 

			}
	if ( object$modtype %in% 2){
		# item parameters
		# l1 <- object$loadings.theta
		cat("\nItem Parameters - Promax Rotated Parameters (THETA)\n",
				"Loadings, Constants, Asymptotes and Descriptives\n\n")
		l1 <- object$promax.theta
		l1 <- data.frame( l1 , "final.constant" = object$final.constants , 
					"lower" = object$lower , "upper"=object$upper,
					"N" = diag(object$N.itempair) , "p" = diag(object$pm) )
		l1 <- round(l1,3)
		print( l1 )
		res <- list( "itempars.promax.theta" = l1 )
		#****
		cat("\nItem Parameters - Promax Rotated Parameters (DELTA)\n",
				"Loadings, Constants, Asymptotes and Descriptives\n\n")
		l1 <- object$promax
		l1 <- data.frame( l1 , "thresh" = object$threshold , 
					"lower" = object$lower , "upper"=object$upper,
					"N" = diag(object$N.itempair) , "p" = diag(object$pm) )
		l1 <- round(l1,3)
		print( l1 )		
				}												
				
	if ( ! is.null(logfile) ){ sink() }				
		invisible(res)
	}
#--------------------------------------------------------------------------------------
