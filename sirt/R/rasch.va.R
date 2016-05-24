
#################################################
# estimating the Rasch model with variational approximation
rasch.va <- function( dat , globconv=.001 , maxiter=1000){
	#********************
	# data preparation
	N <- nrow(dat)
	I <- ncol(dat)
	dat2 <- dat
	dat2.resp <- 1 - is.na(dat) 
	dat2[ dat2.resp==0 ] <- 0
	# initial values
	b <- - stats::qlogis( colMeans( dat , na.rm=T ) )
	sig2 <- 1
	mu.i <- rep(0,N)
	mu.i <- .8*stats::qnorm( rowMeans( ( dat + .1 ) / 1.2 ) )
	sigma2.i <- rep( sig2 , N )
	disp <- "...........................................................\n"	
	iter <- 0
	conv <- 1000
	
	#***************
	# begin algorithm
    while( (globconv < conv)  &	( iter < maxiter )	){
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )		
		b0 <- b
		# compute xsi.ij: formula (18)
		xsi.ij <- matrix( b^2 , N , I , byrow=T) -       # (x'beta)^2
					2 * matrix( b , N , I , byrow=T ) * mu.i +  # 2nd term
					mu.i^2 + sigma2.i             
		xsi.ij <- sqrt( xsi.ij )             
		lam.xsi.ij <- base::tanh( xsi.ij / 2 ) / ( 4*xsi.ij )
		# update beta
		t1 <- colSums( 2 * lam.xsi.ij  * dat2.resp )
		b <- colSums( - ( dat2 - .5 ) * dat2.resp + 2 * lam.xsi.ij * mu.i ) / t1
		# update sigma2
		sig2 <- mean( mu.i^2 + sigma2.i )
		# update sigma2.i
		sigma2.i <- 1 / ( 2 * rowSums( lam.xsi.ij ) + 1 / sig2  )
		# update mu.i
		mu.i <- rowSums( ( dat2 - 1/2 + 2*lam.xsi.ij * matrix( b , N , I , byrow=TRUE ) ) 
								* dat2.resp  ) * sigma2.i
		# display progress
		conv <- max( abs(b-b0))
		iter <- iter+1
		cat( paste( "    Maximum b parameter change = " , 
				paste( round(max(abs(b-b0)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    SD Trait = " , 
				paste( round(sqrt(sig2) ,3) , collapse=" " ) , "\n" , sep=""))				
		utils::flush.console()		
    }
	item <- data.frame("item"=colnames(dat) , "b"=b )
	res <- list("sig" = sqrt(sig2) , 
			"item"=item , "xsi.ij"=xsi.ij , "mu.i"=mu.i , "sigma2.i" = sigma2.i )
	return(res)
		}