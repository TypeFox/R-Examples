mice.impute.pmm4 <- function (y, ry, x, donors=3 , noise = 10^5 , 
		ridge = 10^(-5) , ...){
    x <- cbind(1, as.matrix(x))
#    parm <- .norm.draw(y, ry, x, ...)
    parm <- .norm.draw3(y, ry, x, ridge=ridge ,  ...)	
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
	GG <- 1000* max( abs( yhatobs[,1] ) , abs( yhatmis[,1] ))
	dfr <- cbind( 1 , 1:nrow(yhatobs) , yhatobs[,1] , y[ry] )
	dfr0 <- cbind( 0 , 1:nrow(yhatmis) , yhatmis[,1] , NA)	
	dfr <- rbind( dfr , dfr0 )
	dfr <- data.frame(dfr[ order(dfr[,3] ) , ])
	colnames(dfr) <- c("obs" , "index" , "yhat" , "y")
	# add some small noise to create unique entries in matrix d0
	d00 <- abs(diff(dfr$yhat))
	fg1 <- min( d00[ d00 > 0 ] )	
	dfr$yhat <- dfr$yhat + stats::runif( nrow(dfr) , 0 , fg1 / noise )
	dfr <- data.frame(dfr[ order(dfr[,3] ) , ])		
	dfr$sortindex <- seq( 1 , nrow(dfr))
		
	#******************
	# create donors
	dfr1 <- dfr[ dfr$obs == 1 , ]
	N1 <- nrow(dfr1)
	for (dd in 1:donors){
	#	dd <- 1
		dfr1[ , paste("donor.p",dd,sep="") ] <- c( dfr1$yhat[-seq(1,dd)] , 
					rep(NA,dd) )
		dfr1[ , paste("donor.m",dd,sep="") ] <- c( rep(NA,dd) , 
					dfr1$yhat[-c( seq( N1 ,N1+1 - dd) ) ]  )	
			}
	dfr1$yhatdonor <- dfr1$yhat
	#*****
	# look at appropriate donors
	dfr$sortindex.obs1 <- cumsum(  dfr$obs ) 
	dfr$sortindex.obs2 <- dfr$sortindex.obs1 + 1
	dfr$sortindex.obs1[ dfr$sortindex.obs1 == 0 ] <- 1
	dfr$sortindex.obs2[ dfr$sortindex.obs2 > N1 ] <- N1
	dfr2 <- dfr[ dfr$obs == 0 , ]	
	dfr2a <- dfr1[ dfr2$sortindex.obs1 , grep( "donor" , colnames(dfr1) ) ]
	dfr2a <- cbind( dfr2a , dfr1[ dfr2$sortindex.obs2 , grep( "donor" , colnames(dfr1) ) ]	 )
	dfr2 <- cbind( dfr2 , dfr2a )
	gc <- grep( "donor" , colnames(dfr2))
	d1 <- dfr2.abs <- abs( dfr2$yhat - dfr2[ , gc ] )
	dfr2.abs[ is.na( dfr2.abs ) ] <- GG
	dfr2.gc <- dfr2[,gc]
	yhatmin <- wmin <- matrix( 0 , nrow=sum(!ry) , ncol=donors )
	for ( dd in seq( 1 , donors ) ){
#		dd <- 1
		wmin[,dd] <- apply( dfr2.abs , 1 , min )
		# do.call does not give any speed gains
#		wmin[,dd] <- do.call( pmin , as.data.frame( dfr2.abs) )
		I1 <- ( dfr2.abs == wmin[,dd] )
		dfr2.abs <- dfr2.abs + GG * I1
		yhatmin[,dd] <- rowSums( dfr2.gc * I1 , na.rm=T) / rowSums(I1)
				}
	yhatM <- matrix( t(yhatmin) , ncol=1 , byrow=FALSE )
	wmin2 <- match( yhatM , dfr1$yhat )
	wmin2 <- dfr1[ wmin2 , "index" ]
	wmin2 <- matrix( wmin2 , ncol=donors , byrow=T )
	wmin <- wmin2
	samp <- sample( 1:donors , sum(!ry) , replace = TRUE )
	ind <- wmin[,1] * ( samp == 1 )
	if (donors >= 2){
		for ( dd in seq( 2 , donors )){ 
			ind <- ind + wmin[,dd] * (samp == dd )
						}
					}
	# calculate imputations
	imp <- (y[ry])[ ind ]
	imp <- imp[ order(dfr2$index) ] 
#	flush.console()
    return(imp)
	}
