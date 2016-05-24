


#*************************************************************
# jackknife R2noharm
R2noharm.jackknife <- function( object , jackunits = NULL ){
	#############################################
	# INPUT:
	# object 	... R2noharm object
	# jackunits ... units for jackknife
	#					can be a number or a vector of jackknife units
	#********************************
	# modtype permitted?
	modtype <- object$modtype
	if (modtype %in% 1:2 ){ 
			stop("Jackknife not allowed for exploratory factor analysis.\n") 
				}
	# define jackknife units
	cat("NOHARM4\nJackknife Estimation of Standard Errors\n")
	Nobs <- object$Nobs
	if ( length(jackunits) > 1){ 
			jackunits <- match( jackunits , unique( jackunits) )
								}
	if ( is.null( jackunits ) ){ jackunits <- seq(1,Nobs) }
	if ( length(jackunits) == 1){ jackunits <- 1 + floor( seq( 1 , Nobs - .01 ) / ( Nobs/jackunits ) ) }
	# define progress bar
	u.jackunits <- unique( jackunits )
	UJ <- max( u.jackunits)
	cat( UJ , "Jackknife Units are used\n" )	
	prbar <- paste( "|" , paste(rep("-" , 4 ), collapse="") , "PROGRESS" , 
					paste(rep("-" , 4 ), collapse="") , "|" , sep="")
	cat(prbar,"\n")
	#***************
	# initialize matrices
	cat("|")
	j.tanaka <- rep(NA, UJ )
	j.rmsr <- rep(NA, UJ )
	if (modtype %in% 3:4){
		j.chisquare_df <- j.rmsea <- j.rmsr
						}
	if ( modtype == 4){ 
		j.difficulties <- matrix(NA, nrow= object$Nitems , ncol= UJ )
							}
	j.uniquenesses <- j.final.constants <- j.thresholds <- j.discriminations <- matrix(NA, nrow= object$Nitems , ncol= UJ )
	j.pm <- matrix( NA , nrow= nrow(object$pm)*ncol(object$pm) , ncol=UJ)
	if ( modtype == 3){
		j.factor.cor <- matrix( NA , nrow = prod( dim( object$factor.cor) ) , ncol = UJ )
		j.loadings <- matrix( NA , nrow = prod( dim( object$loadings) ) , ncol = UJ )		
		j.loadings.theta <- matrix( NA , nrow = prod( dim( object$loadings.theta) ) , ncol = UJ )		
						}
	# calculate jackunits indices when progress should be displayed
	prbar.display <- floor( seq( 1 , max(u.jackunits)  , len = 16 ) )
	#**************************
	# BEGIN JACKKNIFE
	for ( uu in u.jackunits ){
		weights.uu <- object$weights[ jackunits != uu ]	
		mod.uu <- R2noharm( dat = object$dat[ jackunits != uu , ] , 
					model.type = object$model.type , 
					weights=weights.uu , 
					guesses = object$guesses , noharm.path = object$noharm.path, 
								F.pattern = object$F.pattern  , F.init = object$F.init , 
								P.pattern = object$P.pattern , P.init = object$P.init , digits.pm = 5 , writename = "model_jackknife" ,
								display.fit = 5  , dec = object$dec , display = FALSE  )
		j.tanaka[uu] <- mod.uu$tanaka
		j.rmsr[uu] <- mod.uu$rmsr
		if (modtype %in% 3:4){
			j.rmsea[uu] <- mod.uu$rmsea
			j.chisquare_df[uu] <- mod.uu$chisquare / mod.uu$df
							}
		if ( modtype == 4){ 
			j.difficulties[,uu] <- mod.uu$difficulties	
			j.discriminations[,uu] <- mod.uu$discriminations
								}
		j.thresholds[,uu] <- mod.uu$thresholds	
		j.uniquenesses[,uu] <- mod.uu$uniquenesses
		j.final.constants[,uu] <- mod.uu$final.constants	
		if (modtype == 3){ 
			j.factor.cor[,uu] <- matrix( mod.uu$factor.cor , ncol=1 )
			j.loadings[,uu] <- matrix( mod.uu$loadings , ncol=1 )			
			j.loadings.theta[,uu] <- matrix( mod.uu$loadings.theta , ncol=1 )						
					}
		j.pm[,uu] <- matrix( as.matrix(mod.uu$pm) , ncol= 1 )
		if ( uu %in% prbar.display ){ cat("-") ; flush.console() }
					}
	cat("|\n")
	#***********************
	# collect jackknife parameters
	jackknife.pars <- list( "j.tanaka" = j.tanaka , "j.rmsr" = j.rmsr )
	if (modtype %in% 3:4){
		jackknife.pars$rmsea <- j.rmsea
		jackknife.pars$chisquare_df <- j.chisquare_df
						}
	if (modtype == 4){ 
		jackknife.pars$j.difficulties = j.difficulties
		jackknife.pars$j.discriminations = j.discriminations 
						}
	jackknife.pars$j.pm = j.pm 
	jackknife.pars$j.thresholds = j.thresholds
	if (modtype == 3){
		jackknife.pars$j.factor.cor <- j.factor.cor
		jackknife.pars$j.loadings <- j.loadings
		jackknife.pars$j.loadings.theta <- j.loadings.theta
					}
	#**********************
	# calculate jackknife statistics

	.est.jackknife.bias <- function( jvec , est , UJ ){
		# newly derived estimators
		if ( is.vector(jvec) ){
			Tn <- UJ * est - (UJ - 1) * jvec
			M <- mean(Tn)
#			se <- sqrt( sum( (Tn - M )^2 ) / (UJ-1 ) )
			se <- sqrt( sum( ( jvec - mean(jvec) )^2 ) * (UJ-1 )/UJ )
							}
		if ( is.matrix(jvec) ){
			Tn <- UJ * est - (UJ - 1) * jvec
			M <- rowMeans(Tn)			
#			se <- sqrt( rowSums( (Tn - M )^2 ) / (UJ-1 ) )
			se <- sqrt(  rowSums( ( jvec - rowMeans(jvec) )^2 ) * (UJ-1 )/UJ )
							}					
		   res <- list( "est.jack" = M , "se" = se)		
				}

			
	res <- .est.jackknife.bias( jvec=j.tanaka , est=object$tanaka , UJ )	
	tanaka.stat <- data.frame( "est" = object$tanaka , 
								"est.jack"=res$est.jack , 
								"se" = res$se 
#								"se" = sqrt( ( UJ - 1 ) / UJ * sum( ( j.tanaka - object$tanaka )^2 ) ) 
										)
	rownames(tanaka.stat) <- "tanaka"
			
	res <- .est.jackknife.bias( jvec=j.rmsr , est=object$rmsr , UJ )			
	rmsr.stat <- data.frame( "est" = object$rmsr , 
								"est.jack"=res$est.jack , 
								"se" = res$se )
	rownames(rmsr.stat) <- "rmsr"
	
	if (modtype %in% 3:4){
	
		res <- .est.jackknife.bias( jvec=j.rmsea , est=object$rmsea , UJ )	
		rmsea.stat <- data.frame( "est" = object$rmsea , 
								"est.jack"=res$est.jack , 
								"se" = res$se )
		rownames(rmsea.stat) <- "rmsea"

        res <- .est.jackknife.bias( jvec=j.chisquare_df , est= object$chisquare / object$df , UJ )	
		chisquare_df.stat <- data.frame( "est" = object$chisquare / object$df , 
								"est.jack"=res$est.jack , 
								"se" = res$se )						
		rownames(chisquare_df.stat) <- "chisquare_df"
			
						}
						
						
	if (modtype == 4){ 


		jvec <- j.difficulties
		est <- object$difficulties
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )	
		difficulties.stat <- data.frame( "est" = est , 
								"est.jack"=res$est.jack , 
								"se" = res$se )	

		jvec <- j.discriminations
		est <- object$discriminations
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )											
		discriminations.stat <-  data.frame( "est" = est , 
								"est.jack"=res$est.jack , 
								"se" = res$se )	
					}

		jvec <- j.thresholds
		est <- object$thresholds 
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )					
		thresholds.stat <- data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )									

		jvec <- j.final.constants
		est <- object$final.constants 
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )									
	final.constants.stat <-data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )
									
		jvec <- j.uniquenesses
		est <- object$uniquenesses 
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )													
	uniquenesses.stat <- data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )
	
	if (modtype == 3){

		jvec <- j.factor.cor
		est <- as.vector(matrix( object$factor.cor , ncol=1 ))
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )			
		factor.cor.stat <- data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )									
		a1a <- rep( rownames(object$factor.cor) , each=ncol( object$factor.cor ) )
		a1b <- rep( colnames(object$factor.cor) , nrow( object$factor.cor ) )
		a1 <- paste( a1a , a1b , sep="-" )
		rownames(factor.cor.stat) <- a1
		
		jvec <- j.loadings
		est <- as.vector(matrix( object$loadings , ncol=1 ))
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )							
		loadings.stat <- data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )								
		a1a <- rep( rownames(object$loadings) , ncol( object$loadings ) )
		a1b <- rep( colnames(object$loadings) , each=nrow( object$loadings ) )
		a1 <- paste( a1a , a1b , sep="-" )
		rownames(loadings.stat) <- a1

		jvec <- j.loadings.theta
		est <- as.vector(matrix( object$loadings.theta , ncol=1 ))
		res <- .est.jackknife.bias( jvec= jvec  , est= est , UJ )					
		loadings.theta.stat <- data.frame( "est" = est , "est.jack"=res$est.jack , "se" = res$se )
		a1a <- rep( rownames(object$loadings.theta) , ncol( object$loadings.theta ) )
		a1b <- rep( colnames(object$loadings.theta) , each=nrow( object$loadings.theta ) )
		a1 <- paste( a1a , a1b , sep="-" )
		rownames(loadings.theta.stat) <- a1		
		
								}
	#*********************************************										
	# collect objects with estimates and standard errors
	se.pars <- list( "tanaka.stat" = tanaka.stat , "rmsr.stat" = rmsr.stat ) 
	if (modtype %in% 3:4){
		se.pars$rmsea.stat <- rmsea.stat
		se.pars$chisquare_df.stat <- chisquare_df.stat
				}
	if (modtype == 4){
		se.pars$difficulties.stat = difficulties.stat 
		se.pars$discriminations.stat = discriminations.stat 
						}
	se.pars$thresholds.stat = thresholds.stat
	se.pars$final.constants.stat = final.constants.stat
	se.pars$uniquenesses.stat = uniquenesses.stat
	if (modtype == 3){
		se.pars$factor.cor.stat <- factor.cor.stat
		se.pars$loadings.stat <- loadings.stat
		se.pars$loadings.theta.stat <- loadings.theta.stat
					}
	# create parameter table				
	res <- se.pars
	LL <- length(res)
	dfr <- NULL
	for (ll in 1:LL ){
		dfr.ll <- res[[ll]]
		dfr.ll <- data.frame( "type" = names(res)[ll] , "parameter" = rownames(dfr.ll) , dfr.ll )
		dfr <- rbind( dfr , dfr.ll )
			}
	rownames(dfr) <- NULL
															
	# output
	res <- list( "partable" = dfr , "se.pars" = se.pars , "jackknife.pars" = jackknife.pars ,
			"u.jackunits"=u.jackunits )
	class(res) <- "R2noharm.jackknife"
	return(res)
	}
#################################################################
	
