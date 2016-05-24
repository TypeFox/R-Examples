 


#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Routine for plausible value imputation
plausible.value.imputation.raschtype <- function( data=NULL , 
		f.yi.qk=NULL , X , 
        Z = NULL, beta0=rep(0,ncol(X)) , sig0=1 ,
        b = rep(1,ncol(X)) ,  a = rep(1 , length(b) ) , c = rep(0 , length(b) ),
        d = 1+0*b , alpha1 =0 , alpha2 =0 ,
        theta.list=seq(-5,5,len=50) , cluster = NULL , 
        iter , burnin , nplausible=1 , printprogress = TRUE ){
    #........................................................................
    # INPUT:                                                            
    # data  ... matrix of dichotomous item responses
    # X     ... matrix of covariates (background variables)
    # Z     ... matrix of covariates allowing for heteroscedasticity
    # beta0 ... initial beta regression coefficients
    # sig0  ... initial residual standard deviation
    # b ... item difficulties
    # a ... item discrimination
    # c ... guessing parameter
    # theta.list ... grid of theta values for posterior density evaluation
    # cluster ... indicates cluster
    # nplausible    ...  number of plausible values
    # iter  ... number of iterations
    # burnin ... number of burn-in iterations    
    #...........................................................................
    # indexes for PV estimation
	if ( ! is.null(data) ){
		data <- as.matrix(data)
		if (is.null(colnames(data)) ){ colnames(data) <- paste( "item" , seq(1 , ncol(data)) , sep="") }
		dfrout <- data.frame( "item" = colnames(data) , "b" = b , "a" = a , "c" = c , "d" = d )  
		cat("\nIRT Plausible Value Imputation - Rasch type model\n")
	    nd <- nrow(data)
						} else {
			nd <- nrow(f.yi.qk) 
					}
    pv.indexes <- sort( unique( round( sort( seq( iter , burnin +1 ,
					len = nplausible ) ) )))
    coefs <- matrix( 0 , nrow=iter , ncol=ncol(X) + ( ! is.null(cluster) ) )
    pvdraws <- matrix( 0 , nrow= length(pv.indexes) , ncol= nd ) 
    # matrices
    X <- as.matrix(X)
    if ( is.null(Z) ){ Z <- matrix( 1 , nrow= nrow(X) , ncol=1) }
	coefsZ <- matrix( 0 , nrow=iter , ncol=ncol(Z))	

	#***************
	# Data preparation
	#***
	# data preparation
	if ( ! is.null(data)){
		y <- data
		y[ is.na(y)] <- 0
		dat2.resp <- 1-is.na(data)
		I <- ncol(y)	# number of items
		dat1 <- data.frame( 1+0*y[,1] , 1 )
					}
	LT <- length( theta.list )
	pi.k <- stats::dnorm( theta.list )
	theta.kM <- matrix( theta.list , nrow= nrow(X) ,
					ncol=length(theta.list)  , byrow=TRUE )
	theta.kM2 <- theta.kM^2
	n <- nrow(X)
	# calculate likelihood if not provided
	if ( is.null( f.yi.qk) ){
		f.yi.qk <- .e.step.raschtype( dat1=dat1 , dat2=y , dat2.resp=dat2.resp , 
			theta.k=theta.list , pi.k=pi.k , I=I , 
			n=nrow(y) , b=b , fixed.a=a , fixed.c=c ,  fixed.d=d , 
			alpha1=alpha1 , alpha2=alpha2 , group=NULL , pseudoll=0 , f.qk.yi=NULL )$f.yi.qk
					}
	post <- f.yi.qk / rowSums( f.yi.qk )
	
    ################################
	# begin iterations
    for (ii in 1:iter){ 
    #   ii <- 1
        # draw plausible values
        if ( ii == 1 ){ X1 <- X }
		########################
		# draw plausible values
		pv1 <- .pv.draw( f.yi.qk , X , beta0 , sig0 , 
                    theta.list , pvdraw = 1 )			
        # calculate adjusted group mean in case of clustering
        if ( !  is.null( cluster ) ){ 
          X1 <- cbind( X , .adj.groupmean( variable = pv1$plausible.value[,1] , 
								cluster ) )
                                } else { X1 <- X }
        # sample latent regression model / draw regression parameters
        s2 <- .sampling.latent.regression.raschtype( pv = pv1$plausible.value[,1] ,
					X = X1 , Z = as.matrix(Z) )
        beta0 <- s2$samp.beta   # sampling of regression coefficients
        sig0 <- s2$fitted.sigma      
        coefs[ii,] <- beta0		
		coefsZ[ii,] <- s2$samp.gamma

        if( ii %in% pv.indexes  ){
            pvdraws[ which( pv.indexes == ii) , ] <- pv1$plausible.value[,1]
                                }
        if ( printprogress ){ 
                if (ii == 1 ){ cat("\n Iteration ") }
                cat(paste(ii , ".",sep="")); 
				utils::flush.console() 
                if ( ( ii %% 10 == 0 ) | ( ii == iter  ) ){ cat("\n     ") }
                        }
                }
	cat("\n")				
    # results pv draw object
    pv1[[  "plausible.value" ]]  <- NULL
    coefs <- coefs[ seq( burnin +1 , iter ) , ]
# "EAP" = EAP , "SE.EAP" = SD.Post	
    res <- list( "coefsX" = coefs , "coefsZ" = coefsZ[seq( burnin +1 , iter ),,drop=FALSE] , 
			"pvdraws" = t(pvdraws) , 
			"posterior" = pv1$posterior.density	 , 
			"EAP"=pv1$EAP , "SE.EAP" = pv1$SE.EAP 	, "pv.indexes" = pv.indexes )
    }
#*******************************************************************************************************
# data=data ; X=X ; weights = rep(1,nrow(data)) ; 
#                                beta.init=rep(0,ncol(X)) ; sigma.init =1 ; max.parchange=.0001 ;
#                                theta.list = seq(-5,5,len=50) ; maxiter = 300 


#***********************************************************************************
# function for drawing plausible values (Raschtype model)
##NS export(plausible.value.draw.raschtype)
.pv.draw <- function( f.yi.qk , X , beta0 , sig0 , 
                    theta.list , pvdraw = 1 ){
	#..................................
    # recode missings
#    y <- data
#    y[ is.na(data) ] <- 1
#    respind <- 1 - is.na(data)
    sig0[ sig0 < 0] <- 0
 #   n <- nrow(y)
    # predicted values from latent regression
    M.Regr <- ( X %*% beta0 )[,1]
    if (length(sig0) > 1){ SD.Regr <- sig0 } else { SD.Regr <- rep( sig0 , n ) }
    # matrix of theta values
	n <- nrow(X)
	l1 <- rep(1,n)
	theta.listM <- thetaM <- matrix( theta.list , nrow=n , ncol=length(theta.list) ,
		byrow=TRUE)
    # compute density resulting from regression
    dens.Regr <- stats::dnorm( thetaM , mean=M.Regr , sd = SD.Regr )
    dens.total <- f.yi.qk * dens.Regr
    dens.total <- dens.total / rowSums( dens.total)	
#   theta.listM <- outer( l1 , theta.list )
    # mean of individual posterior distribution
    EAP <- rowSums( theta.listM * dens.total )
    # SD of posterior distribution
    SD.Post <- sqrt( rowSums( theta.listM^2 * dens.total ) -  EAP^2 )
    # one draw of plausible values
    if (  ! pvdraw ){ pvdraw <- NULL } else { 
            pvdraw <- matrix( stats::rnorm( n*pvdraw , 
				mean = rep(EAP,each=pvdraw) , sd = rep(SD.Post,each=pvdraw) ) , ncol=pvdraw , byrow=T )
                         }
    # results
    res <- list( "theta.grid" = theta.list , "posterior.density" = dens.total , 
		"EAP" = EAP , "SE.EAP" = SD.Post ,
        "plausible.value" = pvdraw , "M.Regr" = M.Regr , "SD.Regr" = SD.Regr )
    return(res)
    }
#***********************************************************************************

#.........................................................................
# sample parameters for latent regression model
.sampling.latent.regression.raschtype <- function( pv , X , 
			Z=rep(1,length(pv)) ){
        # INPUT:
        # pv        ... draw of plausible values
        # X         ... matrix of covariates for latent regression model
        #                 intercept is not automatically included => create vector of ones!!
        # Z         ... matrix of covariates for explaining residual variance
        #.............................................................
        # latent regression model
        mod <- stats::lm( pv ~ 0 + X )
        res <- list( "est.beta" = coef(mod) , "vcov.beta" = vcov(mod) )
        # sample beta parameter
        res$samp.beta <- MASS::mvrnorm( mu = res$est.beta , Sigma = res$vcov.beta )
        # residual standard deviation
        n <- nrow(X) ; p <- ncol(X)
        res$est.sigma <- summary(mod)$sigma
        residuals.mod <- ( stats::resid(mod) )^2   * (n-1) / ( n - p - 1)
        mod1 <- stats::lm( residuals.mod ~ 0 + Z )
#        summary(mod1)
        # sample gamma coefficients for heteroscedasticity
        res$samp.gamma <- MASS::mvrnorm( mu = stats::coef(mod1) , Sigma = stats::vcov(mod1) )
        res$fitted.sigma <- sqrt( stats::fitted(mod1) )
        res$lm.latent.regression <- mod
        res$lm.residualsd <- mod1
        return(res)
        }
#..............................................................................


        #*************************************************************
        # function to calculate adjusted mean: eliminate score on the individual
        # variable <- pv1$plausible.value[,1]
        .adj.groupmean <- function( variable , cluster ){
            a1 <- stats::aggregate( variable ,    list( cluster ) , mean  )
            a2 <- stats::aggregate( 1+0*variable ,    list( cluster ) , sum  )
            ind <- match( cluster , a1[,1] )
            ( a2[ind,2] * a1[ ind , 2] - variable ) / a2[ ind , 2] 
                    }
	#*****************************************
