mice.impute.pmm5 <- function (y, ry, x, donors=3 , noise = 10^5 , 
		ridge = 10^(-5) , ...){
    x <- cbind(1, as.matrix(x))
    parm <- .norm.draw3(y, ry, x, ridge=ridge )   
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
    GG <- 1000* max( abs( yhatobs[,1] ) , abs( yhatmis[,1] ))
    dfr <- cbind( 1 , 1:nrow(yhatobs) , yhatobs[,1] , y[ry] )
    dfr0 <- cbind( 0 , 1:nrow(yhatmis) , yhatmis[,1] , NA)  
    dfr <- rbind( dfr , dfr0 )	
    colnames(dfr) <- c("obs" , "index_obs_miss" , "yhat" , "y") 
    # add some small noise to create unique entries in matrix d0
    d00 <- abs(diff(dfr[,"yhat"]))
    fg1 <- min( d00[ d00 > 0 ] )    
    dfr[,"yhat"] <- dfr[,"yhat"] + stats::runif( nrow(dfr) , 0 , fg1 / noise )
    dfr <- data.frame(dfr[ order(dfr[,3] ) , ])             
    dfr$sortindex <- seq( 1 , nrow(dfr))
	dfr$obsindex_low <- cumsum( dfr$obs )
    ind <- seq( nrow( dfr) , 1 , -1 )
    Ny <- sum( ry)
    N0 <- sum( ! ry )
    c1 <- Ny - cumsum( dfr$obs[ ind ] )   + 1    
#    dfr$obsindex_upp <- c1[ ind ]
    dfr$obsindex_upp <- c1[ ind ]	
	dfr$obsindex_low <- mice::squeeze( dfr$obsindex_low , c(1,Ny)) 
	dfr$obsindex_upp <- mice::squeeze( dfr$obsindex_upp , c(1,Ny)) 	
#    dfr[ dfr$obsindex_low < 1  , "obsindex_low" ] <- 1
#    dfr[ dfr$obsindex_upp > Ny  , "obsindex_upp" ] <- Ny
    dfr0 <- dfr[ dfr$obs == 0 , ]
    dfr1 <- dfr[ dfr$obs == 1 , ]

    # create matrix for sampling
    ydonors <- matrix( NA , nrow=nrow(dfr0) , ncol=2*donors )

	dfr0 <- dfr0[ order(dfr0$index_obs_miss) , ]
    for ( dd in 1:donors){
        ydonors[,dd] <- dfr1[ mice::squeeze( dfr0$obsindex_low - dd + 1 ,c(1,Ny) ) , "y"]
        ydonors[,dd+donors] <- dfr1[ mice::squeeze( dfr0$obsindex_upp + dd - 1 ,c(1,Ny) ) , "y"]    
                        }
    ind.sample <- sample( 1:(2*donors) , N0 , replace = TRUE )
    imp <- ydonors[ cbind( 1:N0 , ind.sample) ]
	# Correction ARb 2013-11-13
#	imp <- imp[ dfr0$index_obs_miss ]
    return(imp)
	}
