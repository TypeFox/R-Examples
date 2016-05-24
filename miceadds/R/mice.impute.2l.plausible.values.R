mice.impute.2l.plausible.values <- function (y, ry, x, type , alpha = NULL  , 
                        alpha.se = 0 ,
                        scale.values = NULL , sig.e.miss = 1000000 , 
						like=NULL , theta=NULL , normal.approx=NULL , 
                        pviter = 15 , imputationWeights = rep(1, length(y)) , 
                        plausible.value.print = TRUE , 
                        pls.facs=NULL , interactions=NULL , quadratics =NULL , ...){  
    #*******
	# old arguments which are now excluded from the function
	itemdiff=NULL ; item.resp = NULL ;
    pvirt.iter = 30 ; pvirt.burnin =10 ; 
	pvirt.printprogress=TRUE ;	
	pvirt.printprogress=FALSE ;	
	
	#********
    # define imputation method
    vname <- get("vname", pos = parent.frame()) # get variable name
    newstate <- get( "newstate" , pos = parent.frame() )  
    pvmethod <- 0
    if ( ! is.null( scale.values[[ vname ]] )){ 
			pvmethod <- 3 
					} 
    if ( ! is.null( like[[vname ]] )){ 
			pvmethod <- 4 			
			     }            
    if (pvmethod == 0){
          if ( ! is.null( alpha[[ vname ]] )){ pvmethod <- 1 }
          if (  is.null( alpha[[vname ]] )){ pvmethod <- 2 }            
                    }  
					
    # define scale type
#    scale.type <- .extract.list.arguments( micearg = scale.type , 
#                           vname = vname , miceargdefault = "parallel" )
	scale.type <- "parallel"
		   
    pls.facs <- .extract.list.arguments( micearg = pls.facs , 
                           vname = vname , miceargdefault = NULL )
    interactions <- .extract.list.arguments( micearg = interactions , 
                           vname = vname , miceargdefault = NULL )
    quadratics <- .extract.list.arguments( micearg = quadratics , 
                           vname = vname , miceargdefault = NULL )
           
    ##############################################################
    # Plausible value imputation according to the Rasch model
	# adapt this to include only the likelihood
	##############################################################
    if (pvmethod == 4){     
         res <- .include.2l.predictors( y=y, x=x , ry=ry , type=type , ... )
         X <- res$X
        #*+*+*+*
        # PLS
        if ( is.null(pls.facs) + is.null(interactions) + is.null(quadratics) < 3 ){
            plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , pls.impMethod = "xplsfacs" , 
                            x = X[,-1] , y = y , ry= rep(TRUE,length(y)) , imputationWeights = imputationWeights , 
                            interactions = interactions, quadratics = quadratics ,  pls.facs = pls.facs ,  ... )$yimp
            X <- plsout[,-1]
                }
        #*+*+* 
         cluster <- res$cluster
         # item response data matrix
		 
		 # extract theta grid and likelihood here!!
         like <- as.matrix( like[[ vname ]] )		 
         theta <- as.matrix( theta[[ vname ]] )		 
		 if ( ! is.null( normal.approx[[ vname ]] ) ){
			normal.approx <- normal.approx 
					} else {
			normal.approx <- TRUE
						}
		
		 X <- X[,-1]	# exclude intercept
		 
		 #-- perform latent regression		 
		 mod0 <- TAM::tam.latreg(like=like, theta=theta, Y = X ,
								control=list( progress=FALSE )  )		 
		 #-- draw plausible values
		 cat("\n")
		 mod1 <- TAM::tam.pv( mod0 , normal.approx=normal.approx ,
		 			nplausible=1 , samp.regr=TRUE 	)
         # pv imputation				  
   	     ximp <- mod1$pv[,2]
	 
                    }
    #############################################
    # Plausible value imputation with known scale scores and standard errors
    if (pvmethod == 3){ 
        M.scale <- scale.values[[ vname ]][[ "M" ]]
        SE.scale <- scale.values[[ vname ]][[ "SE" ]]
        # compute true variance 
        true.var <- var.ytrue <- stats::var( M.scale , na.rm=T)  - mean( (SE.scale[ ! is.na(M.scale) ])^2 , na.rm=T )
        miss <- ( is.na(M.scale) ) | ( is.na(SE.scale ) )
        M.scale[miss] <- Mscale <- mean( M.scale , na.rm=TRUE )
        SE.scale[miss] <- sig.e.miss
        # calculate initial means and variances of posterior distribution
        EAP <- ( SE.scale^(-2)*M.scale + true.var^(-1)*Mscale )/( SE.scale^(-2) + true.var^(-1) )
        Var.EAP <- 1 / ( SE.scale^(-2) + true.var^(-1) )  
		x1 <- x
        # group mean where the actual observation is eliminated
        if ( sum( type == -2 ) > 0 ){
            x1b <- cbind( x[ , type== -2 ] , y )
            gm <- mice.impute.2l.groupmean.elim(y = y , ry = FALSE * ry , x = x1b, type = c(-2,1) )						
            x <- x1 <- cbind( x1 , gm )
			type <- c( type , 1 )
			i2 <- which( type == -2 )
			x <- x[ , -i2]
			type <- type[-i2]
						
                    }       
        # group level predictors
#        res <- .include.2l.predictors( y=y, x=x , ry=ry , type=type , ... )
#        X <- res$X[,-1]
		X <- x
        #*+*+*+*
        # PLS
        if ( is.null(pls.facs) + is.null(interactions) + is.null(quadratics) < 3 ){
			if ( is.null(interactions) ){ 
				interactions <- names(type)[ type == 4 ]
				}			
				plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , 
				      pls.impMethod = "xplsfacs" , 
                      x = X , y = y , ry= rep(TRUE,length(y)) , imputationWeights = imputationWeights , 
                      interactions = interactions, quadratics = quadratics ,  pls.facs = pls.facs ,  ... )$yimp
            X <- plsout[,-1]
                }
				
        #*+*+* 
#        cluster <- res$cluster
        xcov1 <- xcov <- X
		xcov1 <- as.matrix(xcov1)
        # begin iterations for drawing plausible values
        for (iter in 1:pviter){ 
            # draw plausible value for individuals
            y.pv <- stats::rnorm( length(EAP) , mean = EAP , sd = sqrt( Var.EAP) )
            # calculate linear regression
            mod <- stats::lm( y.pv ~ xcov1 ) 
#			 mod <- lm( y.pv ~ xcov1a ) 	
            # draw regression parameters
            v <- stats::vcov(mod)
            beta.star <- stats::coef(mod) + MASS::mvrnorm( 1, mu = rep(0,nrow(v) ) , Sigma = v )
            # calculate residual variance in regression
            sigma2 <- mean( stats::residuals(mod)^2 )
            # fitted regression coefficients
            yfitted <- cbind(1,xcov1) %*% stats::coef(mod)
            # update posterior distribution
            EAP <- ( SE.scale^(-2)*M.scale + sigma2^(-1)*yfitted )/( SE.scale^(-2) + sigma2^(-1) )
            Var.EAP <- 1 / ( SE.scale^(-2) + sigma2^(-1) ) 
            # draw plausible value
            y.pv <- stats::rnorm( length(y) ,  mean = EAP , sd = sqrt( Var.EAP) )
			if ( pvirt.printprogress ){  cat("*");  utils::flush.console()  }
            # add mean plausible value
            if ( sum( type==-2) ){ 
                x1b <- cbind( x[ , type== -2 ] , y.pv )
                gm <- mice.impute.2l.groupmean.elim(y = y , ry = FALSE * ry , x = x1b, type = c(-2,1) )
                xcov1 <- cbind( xcov , gm )
                            }
                    }
    ximp <- as.vector( y.pv )
                }
    #############################################
    # PV imputation scale score according to CTT (parallel measurements)
    #    alpha is known or unknown
    if ( pvmethod %in% c(1,2) ){
        # extract scale values
        if ( sum(type==3) == 0){ 
                cat( "\n",paste( "Items corresponding to scale" , vname , 
                        "must be declared by entries of 3 in the predictor matrix"),"\n")
                               }
        dat.scale <- x[ , type == 3 , drop=FALSE ]
        x1 <- x[ , type %in% c(1,2) ]
        # group level predictors
        res <- .include.2l.predictors( y=y, x=x , ry=ry , type=type , ... )
        x1 <- res$X[,-1]
    #*+*+*+*
    # PLS
    if ( is.null(pls.facs) + is.null(interactions) + is.null(quadratics) < 3 ){
        plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , pls.impMethod = "xplsfacs" , 
                        x = x1 , y = y , ry= rep(TRUE,length(y)) , imputationWeights = imputationWeights , 
                        interactions = interactions, quadratics = quadratics ,  pls.facs = pls.facs ,  ... )$yimp
        x1 <- plsout[,-1]
            }
    #*+*+* 

        cluster <- res$cluster
        # group mean where the actual observation is eliminated
        if ( sum( type == -2 ) > 0 ){
            x1b <- cbind( x[ , type== -2 ] , y )
            gm <- mice.impute.2l.groupmean.elim(y = y , ry = FALSE * ry , x = x1b, type = c(-2,1) )
            x1 <- cbind( x1 , gm )
                    }       
        # compute scale score
        y1 <- rowMeans( dat.scale ) 

        #*******
        # plausible value imputation if alpha is estimated or known
        if (pvmethod  %in% c(1,2) ){
            if (pvmethod == 2){  
                alpha.est <- .cronbach.alpha( dat.scale )
            cirel.type <- "Normal Theory"
#            if (scale.type == "congeneric"){ cirel.type <- "Factor Analytic" }
            cir <- MBESS::ci.reliability( data = dat.scale , type = cirel.type , interval.type = TRUE )
            alpha.est <- cir$Estimated.reliability
            alpha.se <- cir$SE.reliability
                                }
            if (pvmethod == 1){
                     alpha.known <- alpha[[vname]]
                     if ( is.list(  alpha.se  ) ){
                            alpha.se <- alpha.se[[ vname ]] 
                                    }
                    if ( is.null( alpha.se ) ){ alpha.se <- 0 }
                    alpha.est <- alpha.known
                            }   
            # sampling of Cronbach's Alpha
            alpha.samp <- stats::rnorm( 1 , mean = alpha.est , sd = alpha.se ) 
            alpha.samp <- min( .99 , max( .01 , alpha.samp ) )      # restriction of range
            ximp <- draw.pv.ctt( y = y1 , dat.scale = dat.scale , x =x1 , alpha = alpha.samp )        
                        }    
            }
			
			
    # print progress
  if ( plausible.value.print){
    cat("\n",vname , " Plausible value imputation ")
    if (pvmethod == 1){ 
            cat(paste("with known Cronbach's Alpha of",alpha.known , 
			     "and known standard error of" , alpha.se , "\n")  ) 
                    }                        
    if (pvmethod == 2){ 
        cat( paste( "with estimated measurement error variance for" , scale.type , "items\n      "))
        cat( paste("estimated Cronbach's alpha of" , round( alpha.est , 3 ) ) , 
		         "(SE =" , round( alpha.se,3) ,") \n")
                        }
    if (pvmethod %in% c(1,2) ){
        cat( paste("        sampled Cronbach's alpha of" , round( alpha.samp , 3 ) ) , "\n")
            }
    if (pvmethod == 3){ cat("with known scale scores and known measurement" ,
		   "error standard deviations") }                        
    if (pvmethod == 4){ cat("using a provided likelihood") }                        
    cat("\n") ; utils::flush.console()
                                }
    
	# return imputed values	
    return(ximp)	
}
