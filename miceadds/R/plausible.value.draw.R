plausible.value.draw <- function( data , X , beta0 , sig0 , b = b , 
                    a = rep(1 , length(b) ) , c = rep(0 , length(b) ) , 
                    theta.list = seq(-5,5,len=40) , pvdraw = 1 ){
    #........................................................................
    # INPUT:                                                            
    # data      ... matrix of dichotomous item responses
    # X         ... matrix of covariates (background variables)
    # beta0     ... estimated beta coefficients
    # sig0      ... estimated residual standard deviation
    #               This can be a vector in case of heteroscedasticity!
    # b     ...     item difficulties
    # a     ...     item discrimination
    # c     ...     guessing parameter
    # theta.list... grid of theta values for posterior density evaluation
    # pvdraw    ...  should plausible values be drawn? == FALSE or input number of PVs
    #..................................
    # recode missings
    y <- data
    y[ is.na(data) ] <- 1
    respind <- 1 - is.na(data)
    sig0[ sig0 < 0] <- 0
    n <- nrow(y)
    # predicted values from lalent regression
    M.Regr <- ( X %*% beta0 )[,1]
    if (length(sig0) > 1){ SD.Regr <- sig0 } else { SD.Regr <- rep( sig0 , n ) }
    # matrix of theta values
    thetaM <- outer( rep(1,n) , theta.list )
    # compute density resulting from regression
    dens.Regr <- stats::dnorm( thetaM , mean=M.Regr , sd = SD.Regr )
    # conditional distribution of item responses
    dens.Resp <- matrix( 0 , nrow=n , ncol= ncol(dens.Regr) )
    for (tt in seq(1 , length(theta.list)) ){
        # tt <- theta.list[1]
        ptt <- outer( rep(1 , n) , c + (1-c)*stats::plogis( a * ( theta.list[tt] - b ) ) )
        dens.Resp[,tt] <- exp( rowSums( respind *y * log( ptt) + respind*( 1-y) * log( 1-ptt)  ) )
        }
#print( cbind(b,a,c))
    dens.total <- dens.Resp * dens.Regr
    dens.total <- dens.total / rowSums( dens.total)
    theta.listM <- outer( rep(1,n) , theta.list )
    # mean of individual posterior distribution
    EAP <- rowSums( theta.listM * dens.total )
    # SD of posterior distribution
    SD.Post <- sqrt( rowSums( theta.listM^2 * dens.total ) -  EAP^2 )
    # one draw of plausible values
    if (pvdraw == FALSE ){ pvdraw <- NULL } else { 
                pvdraw <- matrix( stats::rnorm( n*pvdraw , mean = rep(EAP,each=pvdraw) , 
						sd = rep(SD.Post,each=pvdraw) ) , ncol=pvdraw , byrow=T )
                         }
    # results
    res <- list( "posterior.density" = dens.total , "EAP" = EAP , "SE.EAP" = SD.Post ,
                "plausible.value" = pvdraw , "M.Regr" = M.Regr , "SD.Regr" = SD.Regr )
    return(res)
    }
