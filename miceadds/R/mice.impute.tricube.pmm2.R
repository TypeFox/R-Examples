mice.impute.tricube.pmm2 <- function (y, ry, x, tricube.pmm.scale= .2 , tricube.boot = FALSE , ...){
	NM <- NULL
    x <- cbind(1, as.matrix(x))	
	# print some informations
    vname <- get("vname", pos = parent.frame()) # get variable name        
	
    t1 <- .extract.list.arguments( micearg = tricube.pmm.scale , 
                           vname = vname , miceargdefault = .2 )
	
    cat( "\n" , paste( vname , "  Imputation Method tricube.pmm with scaling factor" , 
			t1 , "\n"))
    parm <- .norm.draw2(y, ry, x, ...)
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
	ymean <- mean( y[ry] )
	#***
	# bootstrap
	if ( tricube.boot ){
		y1 <- y[ry] ; x1 <- x[ry,]
		B <- length(y1)
		ind <- sample( 1:B , replace = TRUE )
		parm2 <- .norm.draw2(y = y1[ind] , ry = rep(TRUE,B) , x = x1[ind,] , ...)	
		yhatmis <- x[!ry, ] %*% parm2$beta		
				}
	#***
	# R2 calculations
		R2.orig <- 1 - sum( ( y[ry] - yhatobs )^2 ) / sum( ( y[ry] - ymean)^2 )
#		r21 <- 1 - sum( ( y[ry] - yhatobs )^2 ) / sum( ( y[ry] - ymean )^2  )	
		R2.samp <- 1 - sum( ( y[ry] - x[ry, ] %*% parm$beta )^2 ) / sum( ( y[ry] - ymean)^2 ) 
		cat( paste( "  R2 (original data): " , round(R2.orig,4)  , "\n"))
		cat( paste( "  R2 (sampled coefficients): " , round(R2.samp,4)  , "\n"))	
		
#		cat( paste( "  R2 (original data): calc2 " , round( r21,4)  , "\n"))				
				
		if ( tricube.boot ){ 
				R2.boot <- 1 - sum( ( x[ry, ] %*% parm2$beta - y[ry] )^2 ) / sum( ( y[ry] - ymean)^2 ) 		
				cat( paste( "  R2 (Bootstrap): " , round(R2.boot,4)  , "\n"))	
					}
	#***
    # extract scale parameter for tricube pmm
    vname <- get("vname", pos = parent.frame()) 
# print(tricube.pmm.scale)	
#    tricube.pmm.scale <- .extract.list.arguments( micearg = tricube.pmm.scale , 
#                           vname = vname , miceargdefault = .2 )
    utils::flush.console()
    # doing tricube pmm
	# distance function
	dg <- d0 <- d <- abs( outer( yhatmis[,1] , yhatobs[,1] , "-" ) )
	# look for donorset
	DM <- max(d)
	# add some small noise to create unique entries in matrix d0
	d00 <- abs(d0)
	fg1 <- min( d00[ d00 > 0 ] )
	d0 <- d0 + matrix( stats::runif( nrow(d0)*ncol(d0) , 0 , fg1/10000000 ) , ncol=ncol(d0) )
	# DONOR1
	rmin1 <- apply( d0 , 1 , min )	
	d0 <- d0 + DM*( d0 == outer( rmin1 , rep(1,ncol(d0)) ) )
	# DONOR2
	rmin2 <- apply( d0 , 1 , min )
	d0 <- d0 + DM*( d0 == outer( rmin2 , rep(1,ncol(d0)) ) )
	# DONOR3
	rmin3 <- apply( d0 , 1 , min )
	d0 <- d0 + DM*( d0 == outer( rmin3 , rep(1,ncol(d0)) ) )
	rmd <- rowMeans( d )
	xd <- rep(1,ncol(d))
	s.tricube <- outer( t1 * rmd , xd )
	d <- d / s.tricube
	d[ d > 1 ] <- 1
	d <- ( 1 - d^3 )^3
	# redefine donorset
	eps <- .000001 
	d <- d  + eps * ( dg == outer( rmin1 , xd ) )
	d <- d  + eps * ( dg == outer( rmin2 , xd ) )
	d <- d  + eps * ( dg == outer( rmin3 , xd )	)
	prob.x <- d / ( rowSums(d) + .000000000001 )
	probcs <- t( sapply( seq(1,nrow(prob.x)) , FUN = function(ii){ cumsum(prob.x[ii,]) } ) )
	probcs2 <- probcs
	# draw a random number
	RU <- stats::runif( nrow(d)  )
	# write a function which calculates the index when random
	# number exceeded cumulated probability the first time
	matr <- probcs2	
	xy <- outer( RU , rep( 1,ncol(matr)) )
	xy <- matr - xy	
	xy[ xy < 0 ] <- 1.01
	indvec <- apply( xy , 1 , which.min )
	indvec <- ifelse( indvec == 0 , NM , indvec )
	x1 <- ( y[ry] )[ indvec ]
    return(x1)
    }
