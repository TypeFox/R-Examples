latent.regression.em <- function( data , X , weights = rep(1,nrow(data)) , 
                                beta.init=rep(0,ncol(X)) , sigma.init =1 , 
                                b = b , a= rep(1 , length(b)) , c= rep(0 , length(b)) , 
                                max.parchange=.0001 , 
                                theta.list = seq(-5,5,len=50) , maxiter = 300 ){
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
    while( ( parchange > max.parchange ) & ( iter < maxiter ) ){
        # estimate posterior density
        pv1 <- plausible.value.draw( data=data , X=X , beta0 =beta0  , sig0 =sig0 ,
                     b = b , a=a , c = c , theta.list = theta.list , pvdraw= FALSE )
        # estimate latent regression model (linear model)
        mod <- stats::lm( pv1$EAP ~ 0 + X , weights = weights)
        cmod <- stats::coef(mod)
        # Calculation of residual sd
        sigma <- sqrt( mean( weights* ( pv1$SE.EAP^2 + stats::resid(mod)^2 ) ) )
        parchange <- max( abs(sigma - sig0) , abs( cmod - beta0) )
        cat( paste("Iteration " , iter,": max parm. change " , round( parchange , 8 ) ,sep="") , " # Regr. Coeff. " , 
                                        as.vector(cmod) , "\n") ; 
		utils::flush.console()
        # parameter update
        sig0 <- sigma ; beta0 <- cmod ; iter <- iter + 1
                    }
    # standard errors for regression coefficients
    V <- ncol(X)    # number of X variables (predictors)
    scoefs <- matrix( 0 , nrow=V , ncol=5 )
    scoefs <- data.frame(scoefs)
    colnames(scoefs) <- c("est" , "se.simple" , "se" , "t" , "p" )
    rownames(scoefs) <- names(beta0)
    scoefs[,1] <- beta0
    for (vv in 1:V){   
        xvv <- X[,vv]
        scoefs[vv,2] <- sqrt( 1 / sum( weights * xvv^2  / sigma^2 ) )
        h1 <- weights * xvv^2  / sigma^2 * ( 1 - pv1$SE.EAP^2 / sigma^2 )    
        h1[ h1<0] <- 0
        scoefs[vv,3] <- sqrt( 1 / sum( h1 ) )    
        }    
    scoefs$t <- scoefs$est / scoefs$se
    scoefs$p <- 2 * ( 1 - stats::pnorm( abs( scoefs$t ) ) )
    if ( ! is.null( colnames(X) ) ){ rownames(scoefs) <- colnames(X) } # use column names of X
    #********
    # list of results
    res <- list( "iterations" = iter - 1, "maxiter" = maxiter , "max.parchange" = max.parchange ,
                "coef" = beta0 , "summary.coef" = scoefs , "sigma" = sigma )
    return(res)
    }

	
	
	
#.........................................................................
# sample parameters for latent regression model
.sampling.latent.regression <- function( pv , X , Z=rep(1,length(pv)) ){
        # INPUT:
        # pv        ... draw of plausible values
        # X         ... matrix of covariates for latent regression model
        #                 intercept is not automatically included => create vector of ones!!
        # Z         ... matrix of covariates for explaining residual variance
        #.............................................................
        # latent regression model
        mod <- stats::lm( pv ~ 0 + X )
        res <- list( "est.beta" = stats::coef(mod) , "vcov.beta" = stats::vcov(mod) )
        # sample beta parameter
        res$samp.beta <- MASS::mvrnorm( mu = res$est.beta , Sigma = res$vcov.beta )
        # residual standard deviation
        n <- nrow(X)
		p <- ncol(X)
        res$est.sigma <- summary(mod)$sigma
        residuals.mod <- ( stats::resid(mod) ) ^2   * (n-1) / ( n - p - 1)
        mod1 <- stats::lm( residuals.mod ~ 0 + Z )
        summary(mod1)
        # sample gamma coefficients for heteroscedasticity
        samp.gamma <- MASS::mvrnorm( mu = coef(mod1) , Sigma = vcov(mod1) )
        res$fitted.sigma <- sqrt( stats::fitted(mod1) )
        res$lm.latent.regression <- mod
        res$lm.residualsd <- mod1
        return(res)
        }
#..............................................................................

