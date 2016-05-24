#########################################################
# calculation of B matrix
.mml.3pl.computeB <- function( E , gammaslope ){
       dimE <- dim(E)
	   I <- dimE[1]
	   maxK <- dimE[2]
	   D <- dimE[3]
	   B <- array( 0 , dim= c(I , maxK , D ) )
	   for (dd in 1:D){
       for (cc in 1:maxK){
			B[ , cc,dd] <- E[,cc,dd,] %*% gammaslope
						}
					}			
		return(B)
			}
#############################################################
# faster function for computation of item loadings
.mml.3pl.computeB.v2 <- function( Edes , gammaslope , E ){
	B <- .Call("mml_3pl_compute_B_rcpp" ,
			Edes , gammaslope , dim(E) , PACKAGE="TAM")$B
	B <- array( B , dim(E)[1:3] )
	return(B)
		}


###########################################################################
# reduced skillspace estimation
.mml.3pl.skillspace <- function( Ngroup, pi.k , 
			delta.designmatrix , G , delta , delta.fixed ,			
			eps=10^(-10) ){		
	Z <- delta.designmatrix	
	delta0 <- delta
	ND <- nrow(delta)
	covdelta <- list(1:G)
# eps <- 1E-100
	for (gg in 1:G){
		ntheta1 <- Ngroup[gg] * pi.k[,gg]
		ntheta1 <- ntheta1 / sum(ntheta1 )	
# Revalpr("round(ntheta1,5)")		
		lntheta <- log(ntheta1+eps)
# Revalpr("lntheta")
		mod <- stats::lm( lntheta ~ 0 + Z , weights = ntheta1 )
		covbeta <- stats::vcov(mod)		
		beta <- stats::coef(mod)
		if ( ! is.null( delta.fixed ) ){
		# delta.fixed: 1st column: parameter index
		#              2nd column: group index
		#              3rd column: parameter value 
		    ind.gg <- which( delta.fixed[ ,2] == gg )
			if ( length(ind.gg) > 0 ){
				beta[ delta.fixed[ind.gg,1] ] <- delta.fixed[ind.gg,3]
									}
									}
		pi.k[,gg] <- exp( Z %*% beta ) / Ngroup[gg]
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
# Revalpr("round(pi.k[,1],5)")		
		delta[,gg] <- beta
		covdelta[[gg]] <- covbeta
					}
					
					
	res <- list( "pi.k"=pi.k , "delta"=delta , 
			"covdelta" = covdelta )			
			}
##########################################################################
# calculation of expected counts
.mml.3pl.expected.counts <- function( datindw , nitems , maxK , ntheta , hwt){
			# calculate expected counts
			n.ik <- array( 0 , dim=c(nitems , maxK , ntheta ) )
			N.ik <- array( 0 , dim=c( nitems , ntheta ) )	
			for (kk in 1:maxK){   # kk <- 1
				dkk <- datindw[[kk]]
				g1 <- crossprod( dkk , hwt ) 
				n.ik[,kk,] <- g1
				N.ik <- N.ik + g1
						}
			res <- list("n.ik"=n.ik , "N.ik" = N.ik )
			return(res)
			}
#####################################################################

