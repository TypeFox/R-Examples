
###############################################################
truescore.irt <- function( A , B , c=NULL , d =NULL ,
	theta=seq(-3,3,len=21) , error=NULL , pid = NULL ,
	h=.001 ){
	#**********
	if ( is.vector(B) ){ B <- matrix( B , ncol=1 ) }
	if ( is.vector(A) ){ A <- matrix( A , ncol=1 ) }
	nB <- ncol(B)
	maxK <- nB - rowSums( is.na( B ))
	I <- nrow(B)
	if ( is.null(c) ){ c <- rep(0,I ) }
	if ( is.null(d) ){ d <- rep(1,I ) }	
	if ( is.null(pid) ){ pid <- 1:(length(theta)) }
	# true score
	truescore <- .ts.irf( A , B , c , d , theta )
	# calculate standard error of true score
	if ( is.null( error ) ){
		truescore.error <- NULL } else {
		truescore1 <- .ts.irf( A , B , c , d , theta + h )
		truescore.error <- sqrt( ( ( truescore1 - truescore ) /  h )^2 ) * error
				}
	percscore <- truescore / sum( maxK )
	percscore.error <- NULL
	if ( ! is.null( error ) ){
		percscore.error <- truescore.error / sum( maxK )
							}
	# optimization function values
	ind <- intersect( which( !is.na( theta ) ) , which( !is.na(percscore) ) )
	x0 <- theta[ind]
	y0 <- percscore[ind]
	minf <- sum( ifelse( maxK==1 , c , 0 ) ) / sum(maxK) 
	maxf <- sum( ifelse( maxK==1 , d , maxK ) ) / sum(maxK)	
	critfct <- function(x) { 		
	    a <- x[1]
		b <- x[2]
		# define fit function
		sum( ( y0 - minf - (maxf - minf) * stats::plogis( a*x0 + b ) )^2 )
	}		
#    m1 <- glm( y0 ~ x0 , family="binomial" , control=list(maxit=4) )	
#	m1c <- coef(m1)
#	h1 <- optim( c(m1c[2] , -m1c[1] ) , critfct)
	h1 <- stats::optim( c( .5 , 0 ) , critfct )
#	fitval <- c( h1$par[1] , h1$par[2] )
#	names(fitval) <- c("a" , "b" )		
	#*****
	# OUTPUT
	res <- data.frame( "pid"=pid ,"truescore" = truescore )		
	res$truescore.error <- truescore.error 
	res2 <- data.frame( "percscore" = percscore )
	res2$percscore.error <- percscore.error 
    res3 <- data.frame( "lower"=minf , "upper"=maxf , "a" = h1$par[1] , "b" = h1$par[2] )
	res <- cbind( res , res2, res3 )
	return(res)
		}
####################################################
.ts.irf <- function( A , B , c , d , theta ){
    TP <- length(theta)
    truescore <- rep(0, length(theta))
	# extract maximum item categories from B parameters
	nB <- ncol(B)
    maxK <- nB - rowSums( is.na( B ))
    I <- nrow(B)
    scoreM <- matrix( 1:max(maxK) , nrow=TP , ncol=max(maxK) , byrow=TRUE )
    for (ii in 1:I){
        # ii <- 21
        prob.ii <- matrix( NA , TP , maxK[ii] )
        for (kk in seq(1 , maxK[ii] ) ){
            # kk <- 1
            prob.ii[,kk] <- exp( theta * A[ii,kk] + B[ii,kk] )
                    }
        prob.ii <- prob.ii / ( rowSums( prob.ii ) + 1 )
        if ( ( maxK[ii] == 1 ) & ( abs(c[ii])+abs(1-d[ii]) > 0 ) ){
            prob.ii <- c[ii] + ( d[ii] - c[ii] ) * prob.ii
                            }
        truescore <- truescore + rowSums( prob.ii * scoreM[ , seq( 1 , maxK[ii] ) , drop=FALSE] )
                }
    return(truescore)
        }
############################################################