


######################################################################
# calculation of pairwise marginal likelihood
.ll.rasch.pml.probit <- function( b , a , sigma , Q ,eps.corr, itempairs  , IP , eps=10^(-14) ){
	# multidimensional case
	if ( ! is.null(Q) ){ 
		t1 <- a^2 * diag( ( Q %*% sigma %*% t(Q) ) )
				}
	# unidimensional case
	if ( is.null(Q) ){
		t1 <- a^2*sigma^2	
					}	
	xi <-   - b  / sqrt( 1 + t1 )						
    xi1 <- xi[ itempairs[,"item1"] ] 
    xi2 <- xi[ itempairs[,"item2"] ] 
	a1 <- a[ itempairs[,"item1"] ] 	
	a2 <- a[ itempairs[,"item2"] ] 
	if ( ! is.null(Q) ){
		t1 <- a1  * a2 * diag( ( Q[ itempairs[,"item1"] , ] %*% sigma %*% t(Q[  itempairs[,"item2"] , ]) ) )
				}
	if (is.null(Q) ){  # 1-dimensional case
		t1 <- a1*a2*sigma^2	
				}
    cor.Sigma <- ( t1 + eps.corr ) / ( 1 + t1 )
	cor.Sigma[ cor.Sigma > 1 ] <- .99
    a1 <- stats::pnorm( xi )

    itempairs$p1.item1 <- a1[ itempairs$item1 ]
    itempairs$p1.item2 <- a1[ itempairs$item2 ]
	itempairs$p11 <- pbivnorm::pbivnorm( x = xi1 , y = xi2 , rho = cor.Sigma )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
#    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
	for (ii in ind1){
		itempairs[ , ii ] <- ifelse( itempairs[,ii] < eps , eps , itempairs[,ii] )
				}
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    res <-list( "ll" = ll , "itempairs" = itempairs , "b" = b, "sigma" = sigma ,
                "ind1" = ind1 , "ind2" = ind2 )
    return(res)
        }
######################################################################
#              respml <- .update.ll.rasch.pml( respml , b , sigma , itempairs  , IP , eps=10^(-14) )
######################################################################
# update calculation of pairwise marginal likelihood
.update.ll.rasch.pml.probit <- function( respml , b , a , sigma , Q ,eps.corr , 
						itempairs  , IP , eps=10^(-14) ){
#	xi <- - b  / sqrt( 1 + a^2*sigma^2 )
	# multidimensional case
	if ( ! is.null(Q) ){ 
		t1 <- a^2 * diag( ( Q %*% sigma %*% t(Q) ) )
				}
	# unidimensional case
	if ( is.null(Q) ){
		t1 <- a^2*sigma^2	
					}	
	xi <-   - b  / sqrt( 1 + t1 )	
    itempairs <- respml$itempairs
    xi1 <- xi[ itempairs[,"item1"] ] 
    xi2 <- xi[ itempairs[,"item2"] ] 
	a1 <- a[ itempairs[,"item1"] ] 	
	a2 <- a[ itempairs[,"item2"] ] 
	if ( ! is.null(Q) ){
		t1 <- a1  * a2 * diag( ( Q[ itempairs[,"item1"] , ] %*% sigma %*% t(Q[  itempairs[,"item2"] , ]) ) )
				}
	if (is.null(Q) ){  # 1-dimensional case
		t1 <- a1*a2*sigma^2	
				}
    cor.Sigma <- ( t1 + eps.corr ) / ( 1 + t1 ) 	
#    cor.Sigma <- ( a1*a2*sigma^2 + eps.corr ) / ( 1 + a1*a2*sigma^2 ) 
#	    Sigma.ii <- matrix( c(1,0,0,1) , ncol=2 )
#    Sigma.ii[1,2] <- Sigma.ii[2,1] <- cor.Sigma[1] 
    a1 <- stats::pnorm( xi )
    itempairs$p1.item1 <- a1[ itempairs$item1 ]
    itempairs$p1.item2 <- a1[ itempairs$item2 ]
    ind10 <- which( b != respml$b )
    ind10 <- c( which( itempairs$item1 %in% ind10  ) , which( itempairs$item2 %in% ind10  )  )
	itempairs[ind10,"p11"] <- pbivnorm( x = xi1[ind10] , y = xi2[ind10] , 
			rho = cor.Sigma[ind10] )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    res <- list( "ll" = ll , "itempairs" = itempairs , "b" = b, "sigma" = sigma ,
                "ind1" = ind1 , "ind2" = ind2 )
    return(res)
        }
######################################################################




######################################################################
# calculation of pairwise marginal likelihood
.ll.rasch.pml.logit <- function( b , sigma , itempairs  , IP , eps=10^(-14) ,
			p.ki , S.ki ){
	K <- length(p.ki)
	xi <- matrix( 0 , I , K )
	xi1 <- xi2 <- matrix( 0 , IP , K )
	cor.Sigma <- rep(0,K)
	Sigma0 <- matrix( c(1,0,0,1) , ncol=2 )
	Sigma.ii <- as.list( rep(1,K))
	for (kk in 1:K){
		xi[,kk] <-   - b * S.ki[kk] / sqrt( 1 + sigma^2 * S.ki[kk]^2 )
		xi1[,kk] <- xi[ itempairs[,"item1"] , kk ] 
		xi2[,kk] <- xi[ itempairs[,"item2"] , kk ] 
			}
    a1 <- stats::pnorm( xi )
	PKI <- outer( rep(1,I) , p.ki )
	a1 <- rowSums( a1 * PKI )
    itempairs$p1.item1 <- a1[ itempairs$item1 ]
    itempairs$p1.item2 <- a1[ itempairs$item2 ]
	#********
	# calculation of matrix of covariances
	phi.matrix <- matrix( 0 , IP , 9)
	PKI.matrix <- phi.matrix
	for (ii in 1:IP){
		for (kk1 in 1:3){
			for (kk2 in 1:3){ 
			#	kk1 <- 2
			#	kk2 <- 3
				index1 <- 3* ( kk1 - 1 ) + kk2
				index2 <- 3*(kk2 - 1) + kk1
				cor.kk <- S.ki[kk1] * S.ki[kk2] * sigma^2 / sqrt( ( 1 + sigma^2 * S.ki[kk1]^2 ) *
									( 1 + sigma^2 * S.ki[kk2]^2 ) ) 
				phi.matrix[ii,index1] <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper=c(xi1[ii,kk1],xi2[ii,kk2] ),
						mean=c(0,0),
						sigma = matrix( c(1,cor.kk,cor.kk,1) , 2 , 2 ) 
									)
				PKI.matrix[ii,index1] <- p.ki[kk1] * p.ki[kk2]
#				print( paste( kk1 , kk2 , index1 , index2 ))
							}
						}
					}
	itempairs$p11 <- rowSums( phi.matrix * PKI.matrix )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    res <- list( "ll" = ll , "itempairs" = itempairs , "b" = b, "sigma" = sigma ,
                "ind1" = ind1 , "ind2" = ind2 , "phi.matrix" = phi.matrix , 
				"PKI.matrix" = PKI.matrix ,
				"p.ki" = p.ki , "S.ki" = S.ki)
    return(res)
        }
######################################################################
######################################################################
# update calculation of pairwise marginal likelihood
.update.ll.rasch.pml.logit <- function( respml , b , sigma , itempairs  , IP , eps=10^(-14) ){
	# grap objects
#	p.ki <- respml[["p.ki"]]
#	S.ki <- respml[["S.ki"]]
		p.ki <- c( .25220 , .58522 , .16257 )
		S.ki <- c( .90793 , .57778 , .36403 )
#print(respml)
    phi.matrix <- respml$phi.matrix
	PKI.matrix <- respml$PKI.matrix
	itempairs <- respml$itempairs
	#*******************
	K <- length(p.ki)
	xi <- matrix( 0 , I , K )
	xi1 <- xi2 <- matrix( 0 , IP , K )
	cor.Sigma <- rep(0,K)
	Sigma0 <- matrix( c(1,0,0,1) , ncol=2 )
	Sigma.ii <- as.list( rep(1,K))
	for (kk in 1:K){
		xi[,kk] <-   - b * S.ki[kk] / sqrt( 1 + sigma^2 * S.ki[kk]^2 )
		xi1[,kk] <- xi[ itempairs[,"item1"] , kk ] 
		xi2[,kk] <- xi[ itempairs[,"item2"] , kk ] 
			}
    a1 <- stats::pnorm( xi )
	PKI <- outer( rep(1,I) , p.ki )
	a1 <- rowSums( a1 * PKI )
    itempairs$p1.item1 <- a1[ itempairs$item1 ]
    itempairs$p1.item2 <- a1[ itempairs$item2 ]
	#********
	# calculation of matrix of covariances
    ind10 <- which( b != respml$b )
    ind10 <- c( which( itempairs$item1 %in% ind10  ) , which( itempairs$item2 %in% ind10  )  )
	for (ii in ind10 ){
		for (kk1 in 1:3){
			for (kk2 in 1:3){ 
			#	kk1 <- 2
			#	kk2 <- 3
				index1 <- 3* ( kk1 - 1 ) + kk2
				index2 <- 3*(kk2 - 1) + kk1
				cor.kk <- S.ki[kk1] * S.ki[kk2] * sigma^2 / sqrt( ( 1 + sigma^2 * S.ki[kk1]^2 ) *
									( 1 + sigma^2 * S.ki[kk2]^2 ) ) 
				phi.matrix[ii,index1] <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper=c(xi1[ii,kk1],xi2[ii,kk2] ),
						mean=c(0,0),
						sigma = matrix( c(1,cor.kk,cor.kk,1) , 2 , 2 ) 
									)
				PKI.matrix[ii,index1] <- p.ki[kk1] * p.ki[kk2]
#				print( paste( kk1 , kk2 , index1 , index2 ))
							}
						}
					}
	itempairs$p11 <- rowSums( phi.matrix * PKI.matrix )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    res <-list( "ll" = ll , "itempairs" = itempairs , "b" = b, "sigma" = sigma ,
                "ind1" = ind1 , "ind2" = ind2 , "phi.matrix" = phi.matrix , 
				"PKI.matrix" = PKI.matrix ,
				"p.ki" = p.ki , "S.ki" = S.ki)
    return(res)
        }
######################################################################
