



############################################################
# Pairwise estimation with itemclusters
##NS export(rasch.pairwise.itemcluster)
rasch.pairwise.itemcluster <- function( dat , itemcluster = NULL ,
			b.fixed=NULL , conv = .00001 , maxiter = 3000 , 
			progress = TRUE , b.init = NULL, zerosum = FALSE){
	CALL <- match.call()		
    if ( is.null(b.init) ){ 
			b.init <- - stats::qlogis( colMeans( dat , na.rm=TRUE ) ) 
				}
	s1 <- Sys.time()
    I <- ncol(dat)
	dat <- as.matrix(dat)
	dat0 <- dat
	dat[ is.na(dat) ] <- 9
    b <- b.init
	if ( ! is.null(b.fixed) ){
	    b[ b.fixed[,1] ] <- b.fixed[,2] 
		b.fixed <- cbind( b.fixed , exp( b.fixed[,2] ) )
		zerosum <- FALSE
			}		
    # create count tables
    Aij <- t( dat == 0 ) %*% ( dat == 1 )
    Aji <- t( dat == 1 ) %*% ( dat == 0 )
	# set some entries to zero for itemclusters
	clusters <- unique( itemcluster[ itemcluster != 0 ] )
	CC <- length(clusters)
	for (cc in clusters){
		icc <- which( itemcluster == cc )
		Aji[icc,icc] <- Aij[icc,icc] <- 0
				}
    nij <- Aij + Aji
    eps0 <- eps <- exp(  b )
    max.change <- 10
    iter <- 1
    while( max.change > conv ){
        b0 <- b
		eps0 <- eps
#		g1 <- sum( nij[ii,] * ( eps0[ii] + eps0 )^(-1) )	
		m1 <- matrix( eps0 , I , I , byrow=TRUE ) + matrix( eps0 , I , I ) 
		g1 <- rowSums( nij / m1 )
		eps <- rowSums( Aij ) / g1 
        b <-  log(  eps )
		# put item parameter constraints
	    if ( ! is.null(b.fixed) ){
	       eps[ b.fixed[,1] ] <- b.fixed[,3] 
					}				
		if (zerosum){
				b1 <- - log(eps)
				b2 <- b1 - mean(b1)
				eps <- exp( - b2 )	
		  		   }										
        max.change <- max(abs( b - b0 ))
        if ( progress ){
                cat( "PL Iter." , iter , ": max. parm. change = " , 
                        round( max.change , 6 ) , "\n")
                flush.console()
                    } 
        iter <- iter + 1               
                }				
        item <- data.frame( "N" = colSums(1 -is.na(dat0)) , "p" = colMeans( dat0 , na.rm=T ) , 
                        "b" =  log(eps) )
		if ( is.null(itemcluster) ){ itemcluster <- rep(0,I) }
		item$itemcluster <- itemcluster
		s2 <- Sys.time()									
        res <- list( "b" = b , "eps" = eps , "iter" = iter , "conv" = conv , "dat" = dat0 ,"item" = item ,
			"fct" = "rasch.pairwise.itemcluster" , "s1"=s1 , "s2"=s2 , CALL =CALL ) 
        class(res) <- "rasch.pairwise"
        return(res)
       }
#-------------------------------------------------------------------
