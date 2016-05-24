
#------------------------------------------------------------------
# NOHARM exploratory factor analysis                      
R2noharm <- function( dat=NULL , pm=NULL , n=NULL ,model.type , weights=NULL , 
			dimensions = NULL , guesses = NULL , noharm.path , 
                            F.pattern = NULL , F.init = NULL , P.pattern = NULL , P.init = NULL , 
							digits.pm = 4 , writename = NULL ,
                            display.fit = 5  , dec = "." , display = TRUE  ){

    # INPUT:
    # dat ... dataframe
    # model.type ... CFA or EFA
    # dimensions ... dimensions for exploratory factor analysis
    # guesses ... fixed guessing paramters
    # noharm.path ... path for NOHARM console (console must have the name NOHARM87.exe)
    # digits.pm ... digits for calculation of product moment matrix
    # writename ... name of NOHARM Output file
	# dec ... "." or "," for decimal separator 
	# display ... display noharm settings (TRUE or FALSE)
    #............................................................................
    # The matrices F.pattern, F.init, P.pattern and P.init correspond to the 
    # definition in the NOHARM manual
#	noharm.path <- shQuote(noharm.path)
    if (model.type == "CFA" ){ dimensions <- ncol(F.pattern) }
	if ( display ){ 
		cat("Multidimensional Normal Ogive by Harmonic Analysis (NOHARM 4) \n" )
		cat("C. Fraser & R. P. McDonald (2012)\n" )
		cat("For more informations please look at \n ")
#		cat( "	http://people.niagaracollege.ca/cfraser/download/nhCLwindl.html\n\n")
		cat( "  http://noharm.niagararesearch.ca/nh4cldl.html\n\n")
		cat( paste( "Path of NOHARMCL file: \n  ** " , noharm.path , "\n" ) )
		cat( paste( "Path of NOHARM input and output files: \n  ** " , getwd() , "\n\n" ) )		
		if (model.type == "EFA"){ 
				cat(paste( "Exploratory Item Factor Analysis with" , dimensions , "Dimensions" ) , "\n\n" )
				} else {             cat("Confirmatory Item Factor Analysis \n\n" ) }
		flush.console()
				}
	########################################################				
	# data input
	if ( ! is.null(dat) ){
		# allow also for input of moment matrix data				
		dat <- as.matrix(dat)				
		I <- ncol(dat)     # number of items
		n <- nrow(dat)     # number of subjects
		if ( is.null(weights) ){
			weights <- rep(1,n) } else {
			weights <- weights / sum(weights) * n 
				}
		
		# calculate raw product moment matrix
		dat9 <- dat
		dat.resp <- is.na( dat)
		dat9[ dat.resp ] <- 9
		# pairwise product moment matrix
	#    BM <- t( dat9 * (1- dat.resp ) ) %*% ( dat9 * (1- dat.resp  ) )
		 BM <- crossprod( dat9 * (1- dat.resp  )*weights  )
		# number of persons on an item pair
	#    NM <- t((1- dat.resp ) ) %*% (  (1- dat.resp  ) )
		NM <- crossprod(  (1- dat.resp)*weights  )
        BM <- round( BM/NM , digits.pm )		
					}
	inputdat <- TRUE
	########################################################
	if ( ! is.null(pm) ){				
	  if ( is.vector(pm) ){	
		I2 <- length(pm)
		I <- sqrt( 2*I2+1/4 ) - .5		
		BM <- matrix( 0 , I , I )
		colnames(BM) <- paste0("I",1:I)
		vv <- 0
		for (ii in 1:I){
#		ii <- 1
			BM[ ii , 1:ii ] <- pm[ seq(vv+1,vv+ii) ]
			vv <- vv + ii
					}
		BM <-  BM + t(BM)
		diag(BM) <- diag(BM)/2
		     } else { 
				pm <- BM 
				I <- ncol(pm)
					}
		weights <- rep(1,n)
		NM <- matrix(n , I,I)
		dat <- BM[1:2,,drop=FALSE]
		inputdat <- FALSE
						}
    # Exploratory analysis requested?
    EX <- 1 * ( model.type == "EFA" )
    # supply starting values
    IV <- 1 * (model.type == "CFA" )
    
	if ( is.null(guesses) ){
		guesses <- rep(0, I )
					}
				
    # arrange NOHARM Input Files
	s1 <- Sys.time()
    noharm.input <- c( paste( "R2noharm Input file" , s1 ), 
                        paste( I , dimensions , n  , 1 ,  EX , IV , 0 , 0 ,  sep=" ") 
                                        )
    # add guessing parameter
    noharm.input <- c( noharm.input , "   "  , guesses , "  " )

    if (model.type == "CFA"){
        # pattern matrices
        noharm.input <- c( noharm.input , 
            apply( F.pattern , 1 , FUN = function(ll){ paste( ll , collapse= " " ) } ) , " " ,
            sapply( 1:dimensions , FUN = function(ss){ paste( P.pattern[ ss , 1:ss ] , collapse= " " ) }  ) , 
            "   "     )
        # add initial matrices
        if ( is.null(F.init) ){ F.init <- .5*F.pattern }
        if ( is.null(P.init) ){ P.init <- 0.1  + diag( .9 , dimensions) }
        noharm.input <- c( noharm.input , 
                apply( F.init , 1 , FUN = function(ll){ paste( ll , collapse= " " ) } ) , " " ,
                sapply( 1:dimensions , FUN = function(ss){ paste( P.init[ ss , 1:ss ] , collapse= " " ) }  ) 
                        )
         }

        # product moment matrix
        pm <- unlist( sapply( 1:I , FUN = function( ii){  BM[ ii , 1:ii ]  } ) ) 
        LPM <- length(pm)
        h1 <- seq( 1 , LPM , 10 )
        h2 <- c( seq( 10 , LPM , 10 ) )
        if (length(h1) != length(h2) ){ h2 <- c( h2 , LPM ) }
        h3 <- data.frame( h1 , h2 )
        pm <- apply( h3 , 1 , FUN = function(ll){ paste( pm[ ll[1]:ll[2] ] , collapse=" " ) } )


    # add product moment matrix to NOHARM input                      
    noharm.input <- c( noharm.input , " " , pm   )                                  
	if (dec == "," ){ noharm.input <- gsub( "\\." , "," , noharm.input ) }
    # path specifications and starting NOHARM
    current.path <- getwd()
    setwd( noharm.path )
    writeLines( noharm.input , "mymodel.inp" )
#    writeLines( "NoharmCL mymodel.inp mymodel.out 800 0.00001" , "noharm_analysis.bat" )   
    writeLines( "NoharmCL mymodel.inp mymodel.out" , "noharm_analysis.bat" )   
    # working without DOSBox
    system( "noharm_analysis.bat" , show.output.on.console = F  )
    # read NOHARM output   
    noharmout0 <- readLines( "MYMODEL.out" )
	if (dec == "," ){ noharmout0 <- gsub( "," , "\\." , noharmout0 ) }	
	noharmout1 <- noharmout0
    noharmout1 <- c( noharmout1 , rep("",3) , "ENDE" )
    # set current working directory
    setwd( current.path )
    if (!is.null( writename ) ){ 
                writeLines( noharmout0 , paste( writename , ".out" , sep="") ) 
                writeLines( noharm.input , paste( writename , ".inp" , sep="") )                 
                }
    # results for 1-dimensional IRT analysis
    if (dimensions == 1 & model.type == "EFA" ){
			modtype <- 1
	    res <- list( "tanaka" = .noharm.tanaka( noharmout1 ) , 
            "rmsr" = .noharm.rmsr( noharmout1 ) ,
            "N.itempair" = NM ,
            "pm" = BM ,
			"weights" = weights , 
            "guesses" = guesses ,
            "residuals" = .noharm.residuals( noharmout1 , I=I , dat=dat) , 
            "final.constants" = .noharm.itemlevel( noharmout1 , "Final Constants" , I=I , dat=dat) ,
            "thresholds" = .noharm.itemlevel( noharmout1 , "Threshold Values" , I=I , dat=dat) ,
            "uniquenesses" = .noharm.itemlevel( noharmout1 , "Unique Variances" , I=I, dat=dat) , 
            "difficulties" = .noharm.itemlevel( noharmout1 , "Vector B" , I=I, dat=dat) , 
            "discriminations" = .noharm.itemlevel( noharmout1 , "Vector A" , I=I, dat=dat)    		
                )  }
    # results for noharm.efa with more than one dimension
    if (dimensions > 1 & model.type == "EFA" ){	
			modtype <- 2
        res <- list( "tanaka" = .noharm.tanaka( noharmout1 ) , 
                "rmsr" = .noharm.rmsr( noharmout1 ) ,
                "N.itempair" = NM ,
                "pm" = BM ,
                "guesses" = guesses ,
                "residuals" = .noharm.residuals( noharmout1 , I=I , dat=dat) , 
                "final.constants" = .noharm.itemlevel( noharmout1 , "Final Constants",I=I, dat=dat) ,
                "factor.cor" = .noharm.correlations( noharmout1 , "Factor Correlations" ,
								"Promax Rotated", dimensions = dimensions , dat=dat) , 
                "thresholds" = .noharm.itemlevel( noharmout1 , "Threshold Values" , I=I, dat=dat) ,
                "uniquenesses" = .noharm.itemlevel( noharmout1 , "Unique Variances" , I=I, dat=dat) , 
                "unrotated" =  .noharm.loadings( noharmout1  , "Factor Loadings" ,
                                "Varimax Rotated Factor Loadings" ,  dimensions = dimensions , I=I, dat=dat) , 
                "varimax" = .noharm.loadings( noharmout1  , "Varimax Rotated Factor Loadings" ,
                                "Varimax Rotated Coefficients of Theta" ,   dimensions = dimensions , I=I, dat=dat)  , 
                "varimax.theta" =  .noharm.loadings( noharmout1  , "Varimax Rotated Coefficients of Theta" ,
                                    "oblique" ,   dimensions = dimensions , I=I , dat=dat) ,
                "promax" = .noharm.loadings( noharmout1  , "oblique" , 
                                "Factor Correlations" ,   dimensions = dimensions , I=I, dat=dat) ,
                "promax.theta" = .noharm.loadings( noharmout1  , "Promax Rotated Coefficients of Theta" , 
										"ENDE" ,   dimensions = dimensions , I=I, dat=dat) 
                )   
			}
			
    # results for noharm.cfa with more than one dimension
    if (dimensions > 1 & model.type == "CFA" ){
			modtype <- 3			
        res <- list( "tanaka" = .noharm.tanaka( noharmout1 ) , 
                "rmsr" = .noharm.rmsr( noharmout1 ) ,
                "N.itempair" = NM ,
                "pm" = BM ,
                "guesses" = guesses ,
                "residuals" = .noharm.residuals( noharmout1 , I=I , dat=dat) , 
                "final.constants" = .noharm.itemlevel( noharmout1 , "Final Constants",I=I, dat=dat) ,
                "factor.cor" = .noharm.correlations( noharmout1 , "Final Correlations" , "Residual", dimensions = dimensions , dat=dat) , 
                "thresholds" = .noharm.itemlevel( noharmout1 , "Threshold Values" , I=I, dat=dat) ,
                "uniquenesses" = .noharm.itemlevel( noharmout1 , "Unique Variances" , I=I, dat=dat) , 
                "loadings" =  .noharm.loadings( noharmout1  , "Factor Loadings" ,
                                "ENDE" ,  dimensions = dimensions , I=I, dat=dat) , 
                "loadings.theta" = .noharm.loadings( noharmout1  , "Final Coefficients of Theta" ,
                                "Final Correlations of Theta" ,   dimensions = dimensions , I=I, dat=dat)  
					)  

				
            if( !is.null( colnames(F.pattern) ) ){  
					colnames(res$loadings) <- colnames(F.pattern) 
                    colnames(res$loadings.theta) <- colnames(F.pattern)       
                    colnames(res$factor.cor) <- rownames(res$factor.cor) <- colnames(F.pattern)       
                                }                                          
                 }
	# CFA 1 dimension
    if ( ( dimensions == 1 ) & ( model.type == "CFA" ) ){
					modtype <- 4				
			if ( is.null(colnames(F.pattern) ) ){ colnames(F.pattern) <- "F1" }
        res <- list( "tanaka" = .noharm.tanaka( noharmout1 ) , 
                "rmsr" = .noharm.rmsr( noharmout1 ) ,
                "N.itempair" = NM ,
                "pm" = BM ,
                "guesses" = guesses ,
                "residuals" = .noharm.residuals( noharmout1 , I=I , dat=dat) , 
                "final.constants" = .noharm.itemlevel( noharmout1 , "Final Constants",I=I, dat=dat) ,
                "thresholds" = .noharm.itemlevel( noharmout1 , "Threshold Values" , I=I, dat=dat) ,
                "uniquenesses" = .noharm.itemlevel( noharmout1 , "Unique Variances" , I=I, dat=dat) , 
                "loadings.theta" = .noharm.loadings( noharmout1  , "Final Coefficients of Theta" ,
                                "Residual" ,   dimensions = dimensions , I=I, dat=dat)  ,
				"factor.cor" = as.data.frame(matrix(1 ,1 , 1 )) , 
                "difficulties" = .noharm.itemlevel( noharmout1 , "Vector B" , I=I, dat=dat) , 
                "discriminations" = .noharm.itemlevel( noharmout1 , "Vector A" , I=I, dat=dat)  , 
                "loadings" = .noharm.loadings( noharmout1  , "Factor Loadings" ,
                                "LORD" ,   dimensions = dimensions , I=I, dat=dat)  
				#				,
#	            "loadings" =  .noharm.loadings( noharmout1  , "Factor Loadings" ,
#                                "LORD" ,  dimensions = dimensions , I=I, dat=dat)  					
                )   
            if( !is.null( colnames(F.pattern) ) ){  
					colnames(res$loadings) <- colnames(F.pattern) 
					rownames(res$factor.cor) <- colnames(res$factor.cor) <- colnames(F.pattern) 
                    colnames(res$loadings.theta) <- colnames(F.pattern)       
                              } }    
	# collect arguments in output	
	res$model.type <- model.type	
#	res$Nobs <- nrow(dat)
#	res$Nitems <- ncol(dat)
	res$Nobs <- n
	res$Nitems <- I
	res$modtype <- modtype
	res$F.init <- F.init
	res$F.pattern <- F.pattern
	res$P.init <- P.init
	res$P.pattern <- P.pattern
	res$dat <- dat
	res$systime <- s1
	res$guesses <- guesses
	res$noharm.path <- noharm.path
	res$digits.pm <- digits.pm
	res$dec <- dec
	res$display.fit <- display.fit
	if ( modtype %in% c(1,4)){ res$dimensions <- 1 }
	if ( modtype %in% c(2,3)){ res$dimensions <- ncol(res$factor.cor) }	
	#***
	# calculate fit statistic
	if ( modtype %in% 2:4){ 
		RM <- res$residuals
		PV <- diag( res$pm ) 
		g1 <- sqrt( outer( PV * (1-PV) , PV * ( 1-PV ) ) )
		rM <- ( RM / g1  )
		zM <- 0.5 * log( 1 + rM ) - 0.5 * log( 1 - rM ) 
		# chi square
		res$chisquare <- X2 <- ( res$Nobs - 3 ) * sum( zM^2 )
		# calculate number of estimated parameters and degrees of freedom
		I <- res$Nitems
		if (modtype %in% 3:4){
			Nestpars <- I + sum( F.pattern == 1 ) + length( unique( intersect( F.pattern , 2:9999 ) ) )
			Nestpars <- Nestpars + sum( diag(P.pattern) == 1 ) + 
							length( unique( diag(P.pattern)[ diag(P.pattern) > 1] ) )
			g2 <- P.pattern[ lower.tri(P.pattern) ]
			Nestpars <- Nestpars + sum( g2 == 1 ) + length( unique( intersect( g2 , 2:9999 ) ))
			res$Nestpars <- Nestpars
							}
		if ( modtype %in% 2){
			cs <- 1:(dimensions-1)
			res$Nestpars <- Nestpars <- I + I*dimensions - sum(cs)
					}
		res$df <- df <- 0.5*I*(I+1) - Nestpars
		res$chisquare_df  <- res$chisquare / res$df
		# calculate RMSEA
		res$rmsea <- rmsea <- sqrt(max(c( (X2 / res$Nobs ) / df - 1/res$Nobs , 0)))
		# calculate p values
		res$p.chisquare <- 1 - pchisq( res$chisquare , df = res$df )
			}
	# display
	if (display){
		cat( paste( "Tanaka Index=" , round(res$tanaka,display.fit) , sep="") , "\n" )
		cat( paste( "RMSR=" , round(res$rmsr,display.fit) , sep="") , "\n" )
		if ( ! inputdat ){
		   cat("\n**** Note that Jackknife does not work for pm input! **** \n")
				}
			}
	if ( ! inputdat ){
		   res$dat <- NULL
			}			
	res$inputdat <- inputdat
	res$upper <- 1+0*res$guess
	res$lower <- res$guess
	class(res) <- "R2noharm"
    return( res )
    }
#----------------------------------------------------------


