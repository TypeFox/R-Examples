
#############################################
# fit partial credit (or Rasch model)
pcm.fit <- function( b , theta , dat ){

	N <- length(theta)
	if ( is.vector(b) ){
		b <- matrix( b , ncol=1 )
					}						
	K <- ncol(b)
	I <- nrow(b)
	b[ is.na(b) ] <- 999
	dat.ind <-  1 - is.na(dat)
    dat[ is.na(dat) ] <- 0

	# create probability matrix
	rprobs <- array( 0 , dim = c( N , K +1 , I ) )
	for (ii in 1:I){
		rprobs[,,ii ] <- outer( theta , seq(0,K) , "*" )
		rprobs[,-1,ii] <- rprobs[,-1,ii] - outer( rep(1,N) , b[ii,] )
		rprobs[,,ii] <- exp( rprobs[,,ii] )
		rprobs[,,ii] <- rprobs[,,ii] / rowSums( rprobs[,,ii] )
					}
								  
	# expected response
	Eni <- array( 0 , dim= c(N,I) )

	for (ii in 1:I){
		   Eni[,ii] <- rowSums( outer( rep(1,N) , 0:(K) )*rprobs[,,ii]  )
			  }                
	# calculate residuals
	Yni <- dat - Eni
	# calculate variances
	Wni <- array( 0 , dim= c(N,I) )
	for (ii in 1:I){
		   Wni[,ii] <- rowSums( ( outer( rep(1,N) , 0:(K) )  - Eni[,ii] )^2 *rprobs[,,ii]  )
			  }          	
	# calculate kurtosis
	Cni <- array( 0 , dim= c(N,I) )
	for (ii in 1:I){
		   Cni[,ii] <- rowSums( ( outer( rep(1,N) , 0:(K) )  - Eni[,ii] )^4 *rprobs[,,ii]  )
			  }    			  
    # standardized residual
	zni <- Yni / sqrt( Wni )
	#************************************
	# item fit statistics
	N.item <- colSums( dat.ind )	
	#--- Outfit
	outfit <- colSums( zni^2 * dat.ind ) / N.item
	itemfit <- data.frame( "item" = colnames(dat) , 
					"outfit" = outfit )
	qi <- sqrt( colSums( dat.ind * Cni / Wni^2  ) / N.item^2 - 1 / N.item )
	itemfit$outfit.t <- ( itemfit$outfit^(1/3) - 1 ) * ( 3 / qi ) + qi / 3				
    #--- Infit
	itemfit$infit <- colSums( dat.ind * Wni * zni^2 ) / colSums( dat.ind * Wni )
	qi <- sqrt( colSums( dat.ind * ( Cni - Wni^2 )  ) / ( colSums(Wni*dat.ind)  )^2  )
	itemfit$infit.t <- ( itemfit$infit^(1/3) - 1 ) * ( 3 / qi ) + qi / 3
	itemfit0 <- itemfit
	#************************************
	# person fit statistics				
	N.item <- rowSums( dat.ind )	
	#--- Outfit
	outfit <- rowSums( zni^2 * dat.ind ) / N.item
	personfit <- data.frame( "person" = 1:N , "outfit" = outfit )
	qi <- sqrt( rowSums( dat.ind * Cni / Wni^2  ) / N.item^2 - 1 / N.item )
	personfit$outfit.t <- ( personfit$outfit^(1/3) - 1 ) * ( 3 / qi ) + qi / 3				
    #--- Infit
	personfit$infit <- rowSums( dat.ind * Wni * zni^2 ) / rowSums( dat.ind * Wni )
	qi <- sqrt( rowSums( dat.ind * ( Cni - Wni^2 )  ) / ( rowSums(Wni*dat.ind)  )^2  )
	personfit$infit.t <- ( personfit$infit^(1/3) - 1 ) * ( 3 / qi ) + qi / 3
	# output
	res <- list("itemfit"=itemfit , "personfit"=personfit)
	return(res)
			}
#############################################################		