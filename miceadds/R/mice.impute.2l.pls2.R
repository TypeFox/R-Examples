mice.impute.2l.pls2 <- function(y, ry, x , type , pls.facs = NULL , 
                                pls.impMethod = "pmm" , 
                                pls.print.progress = TRUE , 
                                imputationWeights = rep( 1,length(y) ) , 
								pcamaxcols = 1E9 , 
                                tricube.pmm.scale=NULL , min.int.cor = 0 , 
								min.all.cor = 0 , N.largest = 0 ,
								pls.title = NULL , print.dims= TRUE , 
								pls.maxcols=5000 ,	... )
       {
        #...........................................................................#
        # INPUT                                                                     #
        # pls.facs          ... number of factors for PLS regression                #
        # pls.interactions  ... include.interactions                                #
        #                       => type == 4                                        #
        # pls.quadratics    ... include quadratic terms?                            #
        #                       => type == 5                                        #
		# type 				... = 6 : for these variables no interactions will be created	#	
        # pls.impMethod     ... method "norm" or "pmm" or "tricube.pmm"             #
        #                           "xplsfacs" -> return predicted X PLS factors    #
        # pls.print.progress    ... print progress of PLS regression estimation     #
        # imputationWeight  ... vector of weights for imputation                    #
        # min.int.cor       ... minimal correlation for inclusion of interaction    #
        #                           effects                                         #
		# min.all.cor 			... minimal correlation for main effects		    #
		# N.largest			... select N.largest correlations						#
		# pls.title 		... title which is displayed 							#
        #...........................................................................#
        # library(pls)	
# a0 <- Sys.time()
#cat("\n---------------- start variable------\n") 
		time1 <- Sys.time()	
		n <- NULL
        imputationWeights  <- nrow(x) * imputationWeights / sum(imputationWeights)
        vname <- get("vname", pos = parent.frame()) # get variable name        
        
        imp.temp <- get( "newstate" , pos = parent.frame() )
   
#          newstate <- list(it = k,  im = i,  co = j, dep = vname,  meth = theMethod,
#                     log = oldstate$log)
        # extract PLS factors
        pls.facs <- .extract.list.arguments( micearg = pls.facs , 
                           vname = vname , miceargdefault = 20 )
        # extract PLS imputation method
        pls.impMethod <- .extract.list.arguments( micearg = pls.impMethod , 
                            vname = vname , miceargdefault = "pmm" )
        # extract scaling factor for scaling factor in tricube weighted estimation
        tricube.pmm.scale <- .extract.list.arguments( micearg = tricube.pmm.scale , vname , 
                                miceargdefault = .2 )
        # define minimal correlation for interactions
        min.int.cor <- .extract.list.arguments( micearg = min.int.cor , vname , 
                                miceargdefault = 0 ) 

# cat("get all arguments") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
								
        if( pls.print.progress  ){ 
                cat("\n-------------------PLS----------------------\n",vname )
				cat(" 'mice.impute.2l.pls2'")	
			if (print.dims){
				cat("\n.......... Dimensions .............")
				cat("\n dim y   " , length(y) )
				cat("\n dim ry  " , length(ry) , " | sum(!ry)" ,  sum(!ry) )					
				cat("\n dim x   " , dim(x) )										
				cat("\n dim type" , length(type) )
				t1 <- table(type)
				cat("\n table(type)\n"  )					
				print(t1)
				cat("\n\n")
							}
						}					
		# include predictor variables with type != 0
		nt <- names(type)[ type != 0 ]
		nt <- intersect( nt , colnames(x) )
		x10 <- x <- x[ , nt]
		use.ymat <- ( ! is.null( dim(y) ) )
		
		# eliminate too small correlations
#			y1 <- matrix( y , nrow= nrow(x) )
#			c1 <- cor( y1[ry,1] , x[ry,] )
#			vars.elim <- which( abs(c1) < min.cor ) 
#			if ( length(vars.elim) > 0 ){ x <- x[ , - vars.elim ] } 			
		x <- as.matrix(x)	
		
        # standardize x
        if (sd(imputationWeights) > 0){ # with weights
                iW <- outer( imputationWeights , rep(1,ncol(x) ) )
                Mx <- colSums( iW * x )  / colSums( iW )
                SDx <- sqrt( colSums( iW * x^2 ) / colSums(iW) - Mx^2  )
                x0 <- x <- ( x - outer( rep(1,n) , t(Mx) ) ) / outer( rep(1,n) , t(SDx) )
                        } else {  # without weights
#				xx <- x						
#                x0 <- x <- scale( x  )
                x0 <- x <- ma.scale2( x  )		
					# standatdization using ma.scale2 which prevents from
					# producing NaN for variables with a SD of zero
#				x <- as.matrix(x)
#				Mx <- colMeans(xx)
#				Nx <- nrow(x)
#				SDx <- sqrt( ( ( colSums( xx^2 ) - Nx * Mx^2 ) / ( Nx-1) ) + 10^(-10) )							
#				lenx <- ncol(x)
#				ind <- which( SDx < 10^(-9) )				
#				if ( length(ind)>0){
#					lenx <- length(ind)
#					x[,ind] <- x0[,ind] <- ( xx[,ind,drop=FALSE] - 
#								matrix( Mx , nrow=Nx , ncol= lenx , byrow=TRUE ) ) /
#										matrix( SDx , nrow=Nx , ncol= lenx , byrow=TRUE )
#										}
                            }
		#***********************
		# include cluster effect: group mean (eliminating the subject under study)
		if ( sum(type==-2) ){
			x1 <- cbind( y , x[ , which( type==-2) ] )
			type1 <- c( 1,-2)
			ximp <- mice.impute.2l.groupmean.elim( y=y , ry=ry , x = x1 , type=type1 )
			x <- as.matrix( cbind( x , ximp ) )
			colnames(x)[ ncol(x) ] <- "y_aggr"
			x10 <- x0 <- x
			}				
		type <- c( type , 1 )	
		#***************

        N <- ncol(x)
        if( pls.print.progress  ){ 
#                cat("\n-------------------PLS----------------------\n",vname )
                cat(" Imputation: " , imp.temp$im , ", Iteration:", imp.temp$it   )
				if ( ! is.null(pls.title)){ cat("\n ",pls.title) }
                cat( "\n\nImputation using Partial Least Squares Regression\n")
                cat( substring(Sys.time(),1) ,"\n" ) ; flush.console() 
                cat( "\n" , paste( ncol(x10) , "Predictor Variables" , names(y) ) , "\n")  
                cat("Used Variables " , paste(colnames(x10),collapse=" ") , "\n" , sep="" )
#                cat("\n" , "Minimal Absolute Correlation of min.cor =" , min.cor , "\n")				
#                cat( "  Kept" , paste( ncol(x) , "Predictor Variables" , names(y) ) , "\n")  
 #               cat("Used Variables " , paste(colnames(x),collapse=" ") , "\n" , sep="") 			
                                   }
        # extract interactions and quadratic terms
        pls.interactions <- names(type)[ type == 4 ]
        pls.quadratics <- names(type)[ type == 5 ]   

# cat("end quadratics") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
	
        ##############################################
        # create no interactions
        if ( is.null(pls.interactions) ){
                if( pls.print.progress ){   
                        cat("\n" , paste("Created no Interactions" , substring( Sys.time() ,1) ) , "\n") ; flush.console() 
                                    }
                        }  else
        # create some interactions
                        {       
                use.int <- intersect( colnames(x) , pls.interactions  )
                N1 <- length(use.int)
				# standardize x
				cx <- colMeans( x )
				xs <- x - outer( rep(1,nrow(x)) , cx )
                if (N1 > 0){
                    # search for interaction variables in predictorMatrix x?
                    ind.int <- sort( which(  colnames(x) %in% use.int ) )
                    dfr0 <- NULL
                     if( pls.print.progress ){ 
						cat("\nCreate Interactions")
                        cat("\n" , "Minimal Absolute Correlation for Interactions of min.int.cor =" , min.int.cor , "\n\n")
									}								
					N1t <- 0 ; N2t <- 0
					# which interactions should not be created?
					dont.int <- which( colnames(x) %in% ( names(type)[ type == 6 ] )  )				
					# create design matrix
					cols <- setdiff( seq( 1,ncol(x) ) , dont.int )
					dfr <- cbind( rep( ind.int , each=length(cols) ),
								rep(cols , length(ind.int) ) )
					dfr <- dfr[ dfr[ , 1 ] != dfr[,2] , ]
					ind <- intersect( which(dfr[ , 1] %in% ind.int)  , which(dfr[ , 2] %in% ind.int  ) )
					dfr1 <- dfr[ ind , ]
					dfr <- rbind( dfr[ setdiff( seq(1,nrow(dfr)),ind) , ] , 
							dfr1[ dfr1[,1]< dfr1[,2] , ])
					dfr <- dfr[ order( dfr[,1] ) , ]
									
					# create interactions
					res <- create_interactions( y_=y[ry] , xobs_=as.matrix( x[ry,] ) , 
						    xall_=as.matrix(x) ,index_int_ = as.matrix(dfr), 
						    min_int_cor_= min.int.cor , maxcols_= min(nrow(dfr),pls.maxcols) )					
							
# res$N_interactions is number 							
					r1 <- res$allcorrs
					r1[ is.na( r1[,1] ) , 1] <- 0
					res$allcorrs <- r1
					
					# total number of interactions
				    N1t <- nrow(res$index_int)
					# retained number of interactions
					N2t <- ncol( res$xint )
					hx <- res$xint
					index_int2 <- res$index_int				
					index_int2 <- index_int2[ res$allcorrs[,2] == 1 , , drop=FALSE]
					if ( N2t > 0 ){
						colnames(hx) <- paste0( "X" , index_int2[,1] , "." , index_int2[,2] )
						x <- cbind( x , hx )
							}
					if (N2t == 0 ){
						 res$allcorrs <- 0 * is.na(res$allcorrs )
									}
				    N1t <- rowsum( 1+0*res$allcorrs[,1] , res$index_int[,1] )
					N1h <- rowsum( res$allcorrs[,2] , res$index_int[,1] )
				    if (pls.print.progress){
					    cat(" ")   
						cat(paste( seq( 1 , nrow(N1t)) , 
							colnames(x[,ind.int]) ,  
									"Created" , N1t[,1] , 
									"Interactions | Kept", N1h[,1] , "Interactions " , 
											"\n") )				
									}
																		
						#***************************
                            if( pls.print.progress ){
									cat("\n")
                                    cat(paste("Created" , sum(N1t[,1]) , "Interactions in Total | " , substring( Sys.time() ,1) ) , "\n")
                                    flush.console() 
                                    cat("Interactions with " , paste(use.int,collapse=" ") , "\n" , sep="")
                                    cat("Kept " , N2t , " Interactions in Total \n" , sep="")
                                    cat("  Minimal Absolute Correlation for Interactions of min.int.cor =" , min.int.cor , "\n")
                                    flush.console() 
									}              
                                   }
                   }
        ##############################################
        # create quadratic terms
#        if ( pls.interactions[1] == TRUE ){ pls.interactions <- colnames(x0) ; pls.quadratics <- colnames(x) }
#        if ( pls.quadratics[1] == TRUE ){ pls.quadratics <- colnames(x0) }
        pls.quadratics <- union( pls.quadratics , pls.interactions )
        use.quad <- unique( intersect( colnames(x0) ,  pls.quadratics  )    )
        # exclude variables from constructing quadratic terms if they only possess 2 values
        h1 <- apply( as.matrix(x0[ ,use.quad]) , 2 , FUN = function(tt){ length( table(tt) ) } )
        pls.quadratics <- intersect( pls.quadratics , use.quad[ h1 > 2 ] )
        if (  length( pls.quadratics ) > 0 ){
            use.quad <- unique( intersect( colnames(x0) ,  pls.quadratics  )    )
#            x <- cbind( x , x0[ , use.quad ] * x0[, use.quad ] )
			# use standardized variables for creating quadratic terms
            x <- cbind( x , xs[ , use.quad ] * xs[, use.quad ] )
            colnames(x) <- paste("x" , 1:(ncol(x)) , sep="")
            if( pls.print.progress ){  
                        cat("\n" , paste("Created" , length(use.quad) ,"Quadratic Terms" , 
								substring( Sys.time() ,1) ) , "\n") ; flush.console() 
                        cat("Quadratic terms of " , paste(use.quad,collapse=" ") , "\n" , sep="") 
						flush.console()
                           }
                   }
# cat("end interactions") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					   
		#---------------
		# look for minimal absolute correlations for all variables
			# eliminate correlations
			
			if ( ! use.ymat ){ c1 <- cor( y[ry] , x[ry,] ) }
									else {  
											# look for correlations of all the dummy variables
								c1 <- cor( y[ry,] , x[ry,] ) 
								c1 <- apply( abs(c1) , 2 , mean , na.rm=TRUE )
											}
            elim.ind <- which( abs(c1) < min.all.cor )
			N11 <- ncol(x)		

			Nelim <- length(elim.ind)
			if ( ( N11 - Nelim <= 1 ) & (N11>2) ){
					elim.ind <- elim.ind[ -c(1:2) ]
								}
			if ( length(elim.ind) > 0){ 
					x <- x[ , - elim.ind , drop=FALSE] 
							}
			N12 <- ncol(x)
			if ( pls.print.progress){ 
				cat("\n" , paste( N11 , " Predictor Variables" , sep="") )
				cat("\n" , "Minimal Absolute Correlation of min.all.cor =" , min.all.cor , "\n")				
				cat( "  Kept" , paste( N12 , "Predictor Variables" , names(y) ) , "\n")  
				flush.console()
									}
									
								
		# look for largest correlations
			if ( ! use.ymat ){ 
					c1 <- stats::cor( y[ry] , x[ry,] ) 
							}	else {  
											# look for correlations of all the dummy variables
									c1 <- stats::cor( y[ry,] , x[ry,] ) 
									c1 <- apply( abs(c1) , 2 , mean , na.rm=T )
											}
		#***---***---***---***---***---***---***---***---											
			if (FALSE){		# exclusion 2013-11-13										
				if ( ( N.largest < 1 ) & ( N12 > 100 ) ){ 
							N.largest <- N12 + round(N.largest) 	# subtract some value if N.largest is negative
							N.largest <- round( N.largest) 
									}
				if ( ( N.largest < 1 ) & ( N12 <= 100 ) ){
						N.largest <- N12 
								}		
					 }
		#***---***---***---***---***---***---***---***---											
		if (N.largest>0){  # begin N.largest				 
			dfr1 <- data.frame( "index" = seq( 1 , ncol(x) ) , "abs.cor" = abs(as.vector(c1)) )
			dfr1 <- dfr1[ order( dfr1$abs.cor , decreasing=TRUE) , ]
			x <- x[ , dfr1[ 1:N.largest , "index" ] ]
			# look if some columns do have complete missing entries or SD of zero
			cmna1 <- which( colMeans( is.na(x))  == 1 )
			cmna2 <- which( apply( x , 2, stats::sd ) == 0 )
			cmna <- unique( sort( c( cmna1 , cmna2 ) ) )
		
			if ( length(cmna) > 0 ){ 
					x <- x[ , - cmna ] 
					N.largest <- ncol(x)
						}
			if ( pls.print.progress){ 
				cat("\n" , paste( N12 , " Predictor Variables" , sep="") )
				cat("\n" , "Select Predictors with" , N.largest , "Largest Correlations\n")	
					flush.console()
									}
					} # end N.largest
									
		if ( is.vector(x) ){ x <- cbind( x , x10[,1:2] ) }
		if ( dim(x)[2]==0 ){ 
				x <-  x10[,1:2]  
				x[,1] <- x[,1] + rnorm( dim(x)[1] , sd = 10^(-20) )
				x[,2] <- x[,2] + rnorm( dim(x)[1] , sd = 10^(-20) )
					}
	
		if ( ncol(x) > pcamaxcols ){
		  a0 <- Sys.time()
		  xdims <- min( pcamaxcols , nrow(x)-2 )
			cat("\nDimension reduction with Principal Components Analysis\n")
			if (pcamaxcols > 1){
						cat("Dimension of X:" , ncol(x) , " -> Dimension reduction to " , xdims , "dimensions\n")
								}
			if (pcamaxcols < 1){
						cat("Dimension of X:" , ncol(x) , " -> Dimension reduction to " , 
									100 * pcamaxcols , "% of total variance\n")
								}
								
								
#			flush.console()			
#			xpca <- princomp( x )
#			xpca <- prcomp( x )			

			xpca <- pca.covridge(x=as.matrix(x))
			varexpl <- xpca$sdev^2 
			varexpl <- cumsum( varexpl / sum( varexpl) * 100 )
			xdims <- which( varexpl > 100*pcamaxcols )[1]
			cat( " ->" , xdims , "extracted dimensions\n")
			cat("Explained variance:" , round( varexpl[ xdims] , 2 )  , " % " )
#			x <- xpca$x[ , 1:xdims , drop=FALSE]
			x <- xpca$scores[ , 1:xdims , drop=FALSE]
		  a1 <- Sys.time() ; cat("\nTime needed:" , a1-a0 , "\n")
#		  flush.console()
			}

# cat("end processing") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					   
			
		#...
        x10 <- x    # copy data set of predictors
        ##############################################
        ##############################################
        # calculate partial least squares regression
#        if ( is.null( pls.facs ) ){ pls.facs <- ncol(x) }
        nfac <- min(pls.facs,ncol(x))      
        yobs <- y[ ry ]
		if (use.ymat){ yobs <- y[ry,] }
        # center y obs
        weight.obs <- imputationWeights[ ry ]
        weight.obs <- length( weight.obs ) / sum(weight.obs) * weight.obs
        yobs <- yobs - weighted.mean( yobs , weight.obs )
        xobs <- x[ ry , ]
        # include imputationWeights here and calculate weight.obs
        # in the regression model, only PLS factors of X are used
		
        if( stats::sd(weight.obs) > 0){
                yobs <- weight.obs * yobs
                xobs <- outer( weight.obs , rep(1,ncol(xobs) ) ) * xobs
                                }
        if( pls.print.progress  ){ 
                 cat( "\n" , paste( ncol(xobs) , " Dimensions" , sep="")  ) 
                 cat( "\n" , paste( nfac , " PLS factors are used" , sep="") )
				 flush.console()
                 if ( pls.facs == 0){ cat( "\n" , "All" , ncol(x) , "predictors are used (no PLS dimension reduction)")}
                            cat("\n\n" ) 
                 }
        if (pls.facs > 0){
			if ( ! use.ymat ){
			# is this necessary?
#				dfr <- data.frame("Y" = yobs , xobs )
				VV <- ncol(xobs)
#				colnames(dfr)[-1] <- paste( "X" , 1:VV , sep="")
#				pls.formula  <- as.formula( paste( "Y ~ " , paste( paste( "X" , 
#                                     1:VV , sep="") , collapse= " + " ) , sep="") )
						} else {
					
#					dfr <- data.frame("Y" = I(yobs) , xobs )
					VV <- ncol(xobs)
#					colnames(dfr)[-1] <- paste( "X" , 1:VV , sep="")
#					pls.formula  <- as.formula( paste( "Y ~ " , 
#                                paste( paste( "X" , 1:VV , sep="") , collapse= " + " ) , sep="") )
								}
#			if (is.null(pls.method )){ 
#				mod <- plsr( pls.formula , data = dfr , ncomp=nfac , validation="none"  )
				# include new function	
				
				mod <- kernelpls.fit2( X= as.matrix(xobs) , Y=matrix(yobs,ncol=1) ,ncomp=nfac) 	
# Reval( "print( str(mod))")				
#				mod <- plsr( pls.formula , data = dfr , ncomp=nfac , validation="none"  , method = pls.method )
#							}
            if( pls.print.progress ){  
#				smod <- summary(mod , digits = 6 ) 
				print( round( 100*mod$R2 , 2 ))
					}
#            flush.console()
            dfr2 <- x
#            colnames(dfr2) <- colnames(dfr)[-1]
#            colnames(dfr2) <- colnames(x)
            # predict PLS scores for data set (y,x)
#            pmod <- predict( mod , dfr2 , type="scores")
#		 	pmod <- predict.kernelpls.fit2( mod , X=x )
		 	pmod <- predict.kernelpls.fit2( mod , X= as.matrix(x) )						
		    x <- cbind(1, as.matrix(pmod))
            x11a <- x
            #******
            if( pls.print.progress ){  cat( "\nPLS estimation finished " , substring(Sys.time(),1) ,"\n" ) ; flush.console() }
                    }
        if ( pls.facs == 0){ x <- cbind( 1 , x ) }
# cat("end PLS") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1			
        #********************************
        # estimate linear regression     
        if (pls.impMethod != "xplsfacs" ){   		
        if ( sd( imputationWeights) > 0  ){   # if there exists a real sample weight vector
                        x <- cbind(1, as.matrix(x))
                        xobs <- x[ry,]
                        yobs <- y[ry]
                        weights.obs <- imputationWeights[ ry   ]
                        weights.obs <- length(weights.obs) * weights.obs / sum( weights.obs )
                        parm <- .weighted.norm.draw( yobs = yobs , xobs = xobs , ry = ry , y = y , x = x ,
                                            weights.obs = weights.obs , ... )   
                                         }  else {
#                               parm <- .norm.draw( y , ry ,x )  
                                parm <- .norm.draw3( y , ry ,x )  
                                            }
									}
# cat("norm draw") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					   									
        #................................
        # Do Imputation
        if( pls.print.progress  ){  
                    cat( "\n" , paste( "Imputation Method " , pls.impMethod , sep="") , "\n" ) 
              if( pls.impMethod == "tricube.pmm" ){
                    cat( paste( "  scaling factor" , tricube.pmm.scale , "\n"))
                                            }
                                    }
        #print( parm$coef )      # parm$coef are regression coefficients based on complete(d) data
        #print( parm$beta )      # parm$beta are sampled regression coefficients
		
		#*****
		# imputation method 'norm' (mice package)
        if (pls.impMethod == "norm" ){ 
             x1 <- x[  ! ry, ] %*% parm$beta + rnorm(sum(!ry)) * parm$sigma
                                    }
		#*****
		# imputation method 'pmm' (mice package)
        if (pls.impMethod == "pmm" ){ 
            yhatobs <- x[ry, ] %*% parm$coef
            yhatmis <- x[!ry, ] %*% parm$beta
            x1 <- apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs, y = y[ry], ... )
                                }
		#*****
		# imputation method 'pmm5' (miceadds package)	
        if (pls.impMethod == "pmm5" ){ 		
#			x1 <- mice.impute.pmm6( y=y, ry=ry, x=x , ridge= 1e-04 , ...)	
			x1 <- mice.impute.pmm5( y=y, ry=ry, x=x ,  ...)	
								}
		# imputation method 'pmm5' (miceadds package)	
#        if (pls.impMethod == "pmm6" ){ 		
#			x1 <- mice.impute.pmm6( y=y, ry=ry, x=x , ridge=10^(-4) ...)	
#								}								
								
		#*****
		# imputation method 'tricube.pmm'								
        if (pls.impMethod == "tricube.pmm" ){ 
            yhatobs <- x[ry, ] %*% parm$coef
            yhatmis <- x[!ry, ] %*% parm$beta
            x1 <- apply(as.array(yhatmis), 1, .tricube.pmm.match , yhat = yhatobs, y = y[ry], 
                            tricube.pmm.scale = tricube.pmm.scale , ... )
                                }  
		#*****
		# imputation method 'tricube.pmm2'																
        if (pls.impMethod == "tricube.pmm2" ){ 
			x1 <- mice.impute.tricube.pmm2(y=y, ry=ry, x=x, tricube.pmm.scale= tricube.pmm.scale )
                                }  
		if (pls.impMethod == "xplsfacs" ){        x1 <- x   }
# cat("end imputation") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					   		
		time2 <- Sys.time()
        if( pls.print.progress ){  
                    cat( "\nMissing Data Draws finished " , substring(Sys.time(),1) ,"\n" ) ; flush.console() 
					cat( "Time elapsed:" , print(time2 - time1 ) , "|" , time1 , "|" , time2 )
                    cat( "\n...................PLS....................................\n") ; flush.console() 
                               }
        return(x1)
}





#*********************************************************************************
# extract list argument
.extract.list.arguments <- function( micearg , vname , miceargdefault ){
        # micearg   ... name of mice argument
        # vname     ... variable name
        # miceargdefault    ... default for this variable
        if( is.list(micearg) ){
            if ( ! is.null(micearg[[vname]] ) ){
                            micearg <- micearg[[vname]]
                                    }  else   { micearg <- miceargdefault }                      
                                   }
        if ( is.null(micearg) )  {    micearg <- miceargdefault }       
        return( micearg )
            }
#*****************************************************************************





#------------------------------------------------------------------------
# auxiliary function for PLS imputation
.aux.pls.imputation <- function( newstate , vname , pls.impMethod , x , y , ry ,  
                imputationWeights = rep( 1 , length(y)) , 
                interactions , quadratics , pls.facs ,
                ... ){
    # interactions and quadratic terms
    interactions <- .extract.list.arguments( micearg = interactions , 
                           vname = vname , miceargdefault = NULL )                                
    interactions <- intersect( interactions , colnames(x))
    quadratics <- .extract.list.arguments( micearg = quadratics , 
                           vname = vname , miceargdefault = NULL )                
    quadratics <- setdiff( intersect( quadratics , colnames(x)) , interactions )
	if ( is.vector(x) ){ x <- matrix( x , ncol=1 ) }
	# define variable type
    type <- rep( 1 , ncol(x)  )
    names(type) <- colnames(x)
    pls.use <- FALSE
    if ( length( interactions)>0 ){ type[ interactions ] <- 4 ; pls.use <- TRUE  }
    if ( length( quadratics )>0 ){ type[ quadratics ] <- 5 ; pls.use <- TRUE  }
    if ( pls.use & is.null(pls.facs) ){ pls.facs <- 10000 }  
    #.*.*.*..*.*.*..*.*.*.
    # PLS imputation if specified
    pls.facs <- .extract.list.arguments( micearg = pls.facs , 
                           vname = vname , miceargdefault = NULL )
    if ( ! is.null(pls.facs) ){ 
            yimp <- mice.impute.2l.pls(y=y, ry=ry, x=x , type=type , pls.facs = pls.facs , 
                                pls.impMethod = pls.impMethod , 
                                imputationWeights = imputationWeights , ... )           
                            } else { yimp <- NULL }
    res <- list( yimp = yimp , pls.facs = pls.facs )
    return(res)  
        }
#------------------------------------------------------------------------

