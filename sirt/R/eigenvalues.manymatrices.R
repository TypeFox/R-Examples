
###################################################################
# computation of eigenvalues of many matrices
eigenvalues.manymatrices <- function( Sigma.all , itermax=10 , maxconv=.001,
	inverse=FALSE ){
    D <- sqrt( ncol(Sigma.all) )
    N <- nrow( Sigma.all )
    lambda <- matrix( 0 , nrow=N , ncol=D )
    U <- matrix( 0 , nrow=N , ncol=D^2 )
    Sigma1 <- Sigma.all
    for (dd in 1:D){
        # dd <- 1
        res <- .eigenvalue.manymatrices( Sigma.all=Sigma1 , itermax=itermax , maxconv=maxconv )
        lambda[,dd] <- res$lambda
        U[ , D*(dd-1 )+1:D ] <- res$z
        z <- res$z
        for (zz in 1:D){
            Sigma1[ , D*(zz-1) + 1:D  ] <- Sigma1[ , D*(zz-1) + 1:D  ] - lambda[,dd] * z[,zz] * z
                    }
               }
    res <- list("lambda" = lambda , "U" = U )
    # calculate determinant
    res$logdet <- rowSums( log( lambda ) )
    res$det <- exp( res$logdet )
	# compute inverse matrix if required
	if (inverse){
		Sigma.inv <- 0*Sigma.all
		for (dd in 1:D){
		for (zz in 1:D){
			#zz <- 1
			#dd <- 1
			z <- U[ , D*(dd-1) + 1:D  ]
			Sigma.inv[ , D*(zz-1) + 1:D  ] <- Sigma.inv[ , D*(zz-1) + 1:D  ] + 1 / lambda[,dd] * z[,zz] * z
						}}
		res$Sigma.inv <- Sigma.inv
			}
    return(res)
        }
###############################################################################
# computation of the largest eigenvalue and the first eigenvector
.eigenvalue.manymatrices <- function( Sigma.all , itermax=10 , maxconv=.001 ){
        D <- sqrt( ncol( Sigma.all) )
        N <- nrow(Sigma.all)
        z <- matrix(1,nrow=N , ncol=D )
        z <- z / sqrt( rowSums( z^2 ) )
        z0 <- z
        ii <- 0
        conv <- 2*maxconv
        while ( ( ii < itermax ) & ( conv > maxconv) ){
            for (dd in 1:D){ #        dd <- 1
                z[,dd] <- rowSums( Sigma.all[ , D*(dd-1) + 1:D ] * z0 )
                        }         
            v1 <- sqrt(rowSums( z^2 ))
            z <- z / v1
            res <- .rowmax.matrix( abs(z) )
            ind <- as.matrix( data.frame( 1:N , res ) )
            si <- sign( z[ ind ] / z0[ind] )
            z <- si*z
            lambda <- si*v1
            conv <- max( abs( z - z0 ) )     
            z0 <- z
            ii <- ii+1    
                        }
        res <- list( "lambda" = lambda , "z" = z )
        return(res)
            }
##########################################################################
# rowwise maximum in a matrix
.rowmax.matrix <- function(matr){
    ind <- rep(1 , nrow(matr) )
    v1 <- matr[,1]
    for (dd in 2:( ncol(matr) ) ){
        ind <- ifelse( matr[,dd] > v1 , dd , ind )
        v1 <- ifelse( matr[,dd] > v1 , matr[,dd] , v1 )        
                            }
    return(ind)
        }
########################################################################