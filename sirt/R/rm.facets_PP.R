
#############################################################
# person parameter estimation in partial credit model
rm.facets.mle <- function( data , a , b , theta , WLE=FALSE , 
     maxiter=20 , maxincr=3 , h = .001 , convP = .0001 , maxval=9.99 ,
	 progress=TRUE ){

	theta0 <- theta
	N <- length(theta)
	I <- ncol(data)
	iter <- 1

	#-------- begin algorithm
	while ( ( conv > convP ) & (  iter <= maxiter  ) ){	
		
		theta0 <- theta
		
		ll0  <- rm.calc.ll.theta( data , a , b , theta )
		llP1 <- rm.calc.ll.theta( data , a , b , theta+h )
		llM1 <- rm.calc.ll.theta( data , a , b , theta-h )
		
		if (WLE){
			llP2 <- rm.calc.ll.theta( data , a , b , theta+2*h )
			llM2 <- rm.calc.ll.theta( data , a , b , theta-2*h )		   
				}
		
		ll1 <- (llP1 - llM1) / (2*h )
		ll2 <- (llP1 - 2*ll0 + llM1) / h^2
				
		incr <- - ll1 / ll2
		if (WLE){		
			ll3 <- ( llP2 - 2*llP1 + 2*llM1 - llM2 ) / (2*h^3 )	
			incr <- - ll1 / ll2 - ll3 / ( 2*ll2^2 )
				}
						
		maxincr <- maxincr / 1.05
		incr <- ifelse( abs(incr) > maxincr , sign(incr)*maxincr , incr )
		theta <- theta + incr	
		theta <- ifelse( abs(theta) > maxval , sign(theta)*maxval , theta )
		conv <- max( abs( theta - theta0) )
		if (progress){
			cat("* Iteration" , iter , ":" , "maximum parameter change" ,
				round( conv, 5 ) , "\n") 
			utils::flush.console();
					}	
		iter <- iter + 1					
			}
	#------------- end algorithm
	
	# output
	se <- sqrt( abs( - 1 / ll2 ) )
	se <-   ifelse( abs(theta) == maxval , NA , se )		
	theta <- ifelse( theta == maxval , Inf , theta )
	theta <- ifelse( theta == - maxval , -Inf , theta )
	res <- data.frame( "est" = theta , "se" = se )	
	return(res)	
		}
#############################################################



#####################################################
# calculate item response probabilities
calc.pcm <- function( theta , a , b , ii ){
    K <- ncol(b)
    N <- length(theta)
    matrK <- matrix( 0:K , nrow=N , ncol=K+1 , byrow=TRUE)
    eta <- a[ii] * theta * matrK - matrix( c(0,b[ii,]) , nrow=N , ncol=K+1 , byrow=TRUE)
    eta <- exp(eta)
    probs <- eta / rowSums(eta, na.rm=TRUE)    
    return(probs)
            }
#####################################################

######################################################
rm.calc.ll.theta <- function( data , a , b , theta ){
    N <- length(theta)
    I <- ncol(data)
    ll0 <- rep(0,N)
    for (ii in 1:I){
        probs.ii <- calc.pcm( theta , a , b , ii )
        res <- rm.calc.ll( probs.ii , data , ii )
        ll0 <- ll0 + res
                }
    return(ll0)
            }
######################################################

###############################################
# calculate individual likelihood for item ii
rm.calc.ll <- function( probs , data , ii ){
    N <- nrow(data)
    eps <- 1E-20
    # probs <- log( probs + eps )
	probs <- log(probs)
	m1 <- cbind( 1:N , data[,ii] + 1 )
    h1 <- probs[ m1 ] 
    h1[ is.na(h1) ] <- 0
    return(h1)
        }
################################################		