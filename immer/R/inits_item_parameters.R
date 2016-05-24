
################################################
# inits item parameters
inits_itempars <- function( dat , prior ){
	maxK <- apply( dat , 2 , max , na.rm=TRUE )
	K <- max(maxK)
    I <- ncol(dat)
    b <- matrix( NA , nrow=I , ncol=K)
	b1 <- colMeans( dat , na.rm=TRUE ) / maxK
	for (ii in 1:I){
		# ii <- 1
		b[ii , seq(1,maxK[ii],1) ] <- b1[ii] + seq( -2 , 2 , length= maxK[ii]  )
					}
					
	if ( ! is.null( prior$b$M ) ){					
		b <- prior$b$M
					}
	if ( ! is.null( prior$a$M ) ){					
		a <- prior$a$M
					}					
	a <- rep(1,I)
	res <- list( b=b , maxK = maxK , K = K , a = a , I=I)
	return(res)
				}
#################################################				