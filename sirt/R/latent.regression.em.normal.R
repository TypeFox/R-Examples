

#----------------------------------------------------------------------------
# Latent Regression: EM algorithm
#    Rasch type model: based on individual likelihood
latent.regression.em.normal <- function( y, X , sig.e , 
								weights = rep(1,nrow(X)) , 
                                beta.init=rep(0,ncol(X)) , sigma.init =1 , 
                                max.parchange=.0001 , 
                                maxiter = 300 ,
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
	n <- nrow(X)
	
	EAP <- EAP.like <- y
	SE.EAP <- SE.EAP.like <- sig.e

	#....................................
	# begin iterations
    while( ( parchange > max.parchange ) & ( iter < maxiter ) ){
		# estimate latent regression model (linear model)
        mod <- stats::lm( EAP ~ 0 + X , weights = weights)
        cmod <- stats::coef(mod)
        # Calculation of SD
        sigma <- sqrt( mean( weights* ( SE.EAP^2 + stats::resid(mod)^2 ) ) )
		# fitted values
		fmod <- stats::fitted(mod)
		sigma.res <- stats::sd( stats::resid(mod) )		
		prec <- 1 / SE.EAP.like^2 + 1 / sigma^2
		EAP <- ( 1 / SE.EAP.like^2 * EAP.like  + 1/sigma^2 * fmod )/prec
		SE.EAP  <- 1 / sqrt( prec )
        parchange <- max( abs(sigma - sig0) , abs( cmod - beta0) )
		if (progress){ 
			cat( paste("Iteration " , iter,": max parm. change " , round( parchange , 6 ) ,sep="") ,
					 "\n") ; utils::flush.console()
					 }
        # parameter update
        sig0 <- sigma ; beta0 <- cmod ; iter <- iter + 1
                    }
	############# END ALGORITHM ######################################
	##################################################################
    # standard errors for regression coefficients
    V <- ncol(X)    	# number of X variables (predictors)
	W <- diag(weights)	# diagonal matrix of weights
	# simple covariance matrix
	h1 <- crossprod( X , W ) %*% X	
	h2 <- solve(h1)
	vcov.simple <- sigma^2 * h2
	# covariance matrix which includes measurement error
	error.weights <- 1 - SE.EAP^2 / sigma^2 
	X1 <- X * sqrt(error.weights)
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
	totalvar <- explvar + sigma^2
	rsquared <- explvar / totalvar
	scoefs$beta <- scoefs$est / sqrt( totalvar ) * apply( X , 2 , stats::sd )
	scoefs$fmi <- 1 -  scoefs$se.simple^2 / scoefs$se^2
#	scoefs[,"N.simple" ] <- nrow(data)
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
				"EAP" = EAP , "SE.EAP" = SE.EAP , 
				"explvar" = explvar , "totalvar" = totalvar ,
				"rsquared" = rsquared )
	class(res) <- "latent.regression"				
    return(res)
    }
#-----------------------------------------------------------------------------------------------------

