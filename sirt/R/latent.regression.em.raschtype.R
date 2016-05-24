

#----------------------------------------------------------------------------
# Latent Regression: EM algorithm
#    Rasch type model: based on individual likelihood
latent.regression.em.raschtype <- function( data=NULL , f.yi.qk=NULL , X , 
								weights = rep(1,nrow(X)) , 
                                beta.init=rep(0,ncol(X)) , sigma.init =1 , 
                                b=rep(0,ncol(X)) , a= rep(1 , length(b)) , c= rep(0 , length(b)) , 
                                d = rep(1 , length(b)) , alpha1=0 , alpha2=0 , 
                                max.parchange=.0001 , 
                                theta.list = seq(-5,5,len=20) , maxiter = 300 ,
								progress=TRUE ){
    #.....................................................................................#
    # INPUT:                                                                              #
    # data          ... item response matrix                                              #
    # X             ... covariates (background variables)                                 #
    # weights       ... sample weights                                                    #
    # beta.init     ... initial estimate of beta coefficients                             #
    # sigma.init    ... initial estimate of residual standard deviation                   #
    # b     ... item difficulties                                                 #
    # a ... item discrimination
    # c ... guessing parameter
    # d ... 1 minus slipping parameter
    # alpha1, alpha2 ... alpha parameter in rasch type model                              #
    # max.parchange ... maximum parameter change                                          #
    # theta.list    ... grid of theta values for evaluation of posterior density          #
    # maxiter       ... maximum number of iterations                                      #
    #.....................................................................................#
    X <- as.matrix(X)   # matrix format
    # init parameters
    beta0 <- beta.init
    sig0 <- sigma.init
    # normalize sample weights
    weights <- length(weights) * weights / sum(weights )
    # initialize iteration index and value of parameter change
    iter <- 1 ; parchange <- 1000
	#***
	# data preparation
	if ( ! is.null(data)){
		data <- as.matrix(data)	
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
	#....................................
	# begin iterations
    while( ( parchange > max.parchange ) & ( iter < maxiter ) ){

		EAP <- rowSums( theta.kM * post )
		SE.EAP <- sqrt( rowSums( theta.kM2 * post ) - EAP^2	 )
		# estimate latent regression model (linear model)
        mod <- stats::lm( EAP ~ 0 + X , weights = weights)
        cmod <- stats::coef(mod)
        # Calculation of SD
        sigma <- sqrt( mean( weights* ( SE.EAP^2 + stats::resid(mod)^2 ) ) )
		# fitted values
		fmod <- stats::fitted(mod)
		sigma.res <- stats::sd( stats::resid(mod) )		
		h1 <- theta.kM - matrix( fmod ,  nrow=n , ncol=LT )
		prior <- stats::dnorm( h1 , sd = sigma )
		prior <- prior / rowSums( prior )
		post <- f.yi.qk * prior
		post <- post / rowSums( post )
        parchange <- max( abs(sigma - sig0) , abs( cmod - beta0) )
		if (progress){ 
			cat( paste("Iteration " , iter,": max parm. change " , round( parchange , 6 ) ,sep="") ,
	#				" # Regr. Coeff. " , round(as.vector(cmod),5) , 
					 "\n") ; flush.console()
					 }
        # parameter update
        sig0 <- sigma ; beta0 <- cmod ; iter <- iter + 1
                    }
    # standard errors for regression coefficients
    V <- ncol(X)    	# number of X variables (predictors)
	W <- diag(weights)	# diagonal matrix of weights
	# simple covariance matrix
#	h1 <- t(X) %*% W %*% X
	h1 <- crossprod( X , W ) %*% X	
	h2 <- solve(h1)
	vcov.simple <- sigma^2 * h2
	# covariance matrix which includes measurement error
	error.weights <- 1 - SE.EAP^2 / sigma^2 
#	EW <- outer( sqrt( error.weights ) , rep(1,ncol(X)))
#	X1 <- X * EW
	X1 <- X * sqrt(error.weights)
#	h1 <- t(X1) %*% W %*% X1
	h1 <- crossprod( X1 , W ) %*% X1	
	h2 <- solve(h1)
	vcov.latent <- sigma^2 * h2
    scoefs <- matrix( 0 , nrow=V , ncol=9 )
    scoefs <- data.frame(scoefs)
    colnames(scoefs) <- c("est" , "se.simple" , "se" , "t" , "p" , 
							"beta" , "fmi" , "N.simple" , "pseudoN.latent" )
    rownames(scoefs) <- names(beta0)
    scoefs[,1] <- beta0
	scoefs[, "se.simple"] <- sqrt( diag( vcov.simple) )
	scoefs[, "se"] <- sqrt( diag( vcov.latent) )
	# explained variance
	explvar <- stats::var( fmod )
#	totalvar <- explvar + sigma^2
	totalvar <- explvar + sigma^2
	rsquared <- explvar / totalvar
	scoefs$beta <- scoefs$est / sqrt( totalvar ) * apply( X , 2 , stats::sd )
	scoefs$fmi <- 1 -  scoefs$se.simple^2 / scoefs$se^2
	scoefs[,"N.simple" ] <- nrow(X)
	scoefs[,"pseudoN.latent"] <- sum( error.weights )  
    scoefs$t <- scoefs$est / scoefs$se
    scoefs$p <- 2 * ( 1 - stats::pnorm( abs( scoefs$t ) ) )
    if ( ! is.null( colnames(X) ) ){ rownames(scoefs) <- colnames(X) } # use column names of X
	if (progress){
		cat("\nRegression Parameters\n\n")
		.prnum(scoefs,4)		# print results
		cat( paste( "\nResidual Variance  =" , round( sigma^2 , 4 ) ) , "\n" )
		cat( paste( "Explained Variance =" , round( explvar , 4 ) ) , "\n" )
		cat( paste( "Total Variance     =" , round( totalvar , 4 ) ) , "\n" )
		cat( paste( "	        R2 =" , round( rsquared , 4 ) ) , "\n" )	
				}
    #********
    # list of results
    res <- list( "iterations" = iter - 1, "maxiter" = maxiter , 
				"max.parchange" = max.parchange ,
                "coef" = beta0 , "summary.coef" = scoefs , "sigma" = sigma ,
				"vcov.simple" = vcov.simple , "vcov.latent" = vcov.latent , 
				"post" = post , "EAP" = EAP , "SE.EAP" = SE.EAP , 
				"explvar" = explvar , "totalvar" = totalvar ,
				"rsquared" = rsquared )
	class(res) <- "latent.regression"
    return(res)
    }
#-----------------------------------------------------------------------------------------------------

