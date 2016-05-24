

###################################################
# fit logistic model
fit.gradedresponse <- function( freq.categories , SC , I , K , 
	conv=.0001 , maxit=100 , progress=TRUE ){
	#*************************  
	if (progress){ cat("\n*******Graded Response Model***********\n") }	
	theta <- stats::qlogis( seq( .5 , SC-1 , len=SC ) / SC )
	# item parameters
	b <- rep(0,I)
	b.cat <- seq(1.5 , -1.5 , len=K)
	#***********************************
	# begin algorithm
    numdiff.parm <- .001
	max.increment <- 1
	increment.factor <- 1.01
	deviation <- 1000
	iter <- 0
	########
	# algorithm
	while( ( iter < maxit) & ( deviation > conv ) ){   
		b0 <- b
		b.cat0 <- b.cat
		theta0 <- theta
		max.increment <- max.increment * increment.factor^(-iter)
		# update theta
		res <- .update.theta.grm( theta , b , b.cat , freq.categories ,
			numdiff.parm , max.increment=max.increment)							
		theta <- res$theta
		ll <- res$ll
		prob <- res$prob		
		# update b
		res <- .update.b.grm( theta , b , b.cat , freq.categories ,
			numdiff.parm , max.increment)			
		b <- res$b
        b <- b - mean(b)
		# update b.cat
		res <- .update.bcat.grm( theta , b , b.cat , freq.categories ,
			numdiff.parm , max.increment)	
		b.cat <- res$b.cat	
		b.cat <- b.cat - mean(b.cat )
	    deviation <- max( abs( c( theta-theta0 , b - b0 , b.cat - b.cat0 )) )
		if (progress){ 
				cat( "Iteration" , iter , "- Deviation =" ,  round( deviation , 6 ) , "\n")
				flush.console()
					}	
        iter <- iter + 1					
			}  # end algorithm
	#*************************
	llcase.grm <- ll / mean( rowSums(colSums( freq.categories )))
	# output
	res <- list( "item.sc"= b  , "cat.sc"=b.cat ,
				"person.sc" = theta , "ll"=ll  , "llcase" = llcase.grm , 
				"prob"=prob )
    return(res)
	}
####################################################################