plausible.value.imputation <- function( data , X , Z = NULL, beta0=rep(0,ncol(X)) , sig0=1 ,
                b = b ,  a = rep(1 , length(b) ) , c = rep(0 , length(b) ),
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
    dfrout <- data.frame( "item" = colnames(data) , "b" = b , "a" = a , "c" = c )  
    cat("\nIRT Plausible Value Imputation 3 PL with item parameters\n")
    print(dfrout)
    pv.indexes <- sort( unique( round( sort( seq( iter , burnin +1 , len = nplausible ) ) )))
    coefs <- matrix( 0 , nrow=iter , ncol=ncol(X) + ( ! is.null(cluster) ) )
    pvdraws <- matrix( 0 , nrow= length(pv.indexes) , ncol=nrow(data) ) 
    # matrices
    X <- as.matrix(X)
    if ( is.null(Z) ){ Z <- matrix( 1 , nrow= nrow(X) , ncol=1) }
        #*************************************************************
        # function to calculate adjusted mean: eliminate score on the individual
        # variable <- pv1$plausible.value[,1]
        .adj.groupmean <- function( variable , cluster ){
            a1 <- stats::aggregate( variable ,    list( cluster ) , mean  )
            a2 <- stats::aggregate( 1+0*variable ,    list( cluster ) , sum  )
            ind <- match( cluster , a1[,1] )
            ( a2[ind,2] * a1[ ind , 2] - variable ) / a2[ ind , 2] 
                    }
        #*************************************************************    
    for (ii in 1:iter){ 
    #   ii <- 1
        # draw plausible values
        if ( ii == 1 ){ X1 <- X }
        pv1 <- plausible.value.draw( data=data , X=X1 , beta0 = beta0  , sig0 =sig0 , 
                        b = b , a = a , c=c , theta.list = theta.list , pvdraw=1 )
        # calculate adjusted group mean in case of clustering
        if ( !  is.null( cluster ) ){ 
                       X1 <- cbind( X ,  .adj.groupmean( variable = pv1$plausible.value[,1] , cluster ) )
                                } else { X1 <- X }
        # sample latent regression model / draw regression parameters
        s2 <- .sampling.latent.regression( pv = pv1$plausible.value[,1] , X = X1 , Z = as.matrix(Z) )
        beta0 <- s2$samp.beta   # sampling of regression coefficients
        sig0 <- s2$fitted.sigma      
        coefs[ii,] <- beta0
        if( ii %in% pv.indexes  ){
            pvdraws[ which( pv.indexes == ii) , ] <- pv1$plausible.value[,1]
                                }
        if ( printprogress ){ 
                if (ii == 1 ){ cat("\n Iteration ") }
                cat(paste(ii , ".",sep="")); utils::flush.console() 
                if ( ( ii %% 10 == 0 ) | ( ii == iter  ) ){ cat("\n     ") }
                        }
                }
    coefs <- coefs[ seq( burnin +1 , iter ) , ]
    res <- list( "coefs" = coefs , "pvdraws" = t(pvdraws) , "pv.indexes" = pv.indexes )
    }
