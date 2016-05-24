

###################################################
# fit logistic model
fit.logistic <- function( freq.correct , wgt , scores , item.p ,
	conv=.0001 , maxit=100 , progress=TRUE ){
	#*************************
    sc.adj <- min( diff(scores) ) / 3
	scores1 <- scores
	freq.correct <- as.matrix(freq.correct)
	RR <- nrow(freq.correct)
	CC <- ncol(freq.correct)
	scores1[ scores == 0 ] <- sc.adj
	scores1[ scores == 1 ] <- 1 - sc.adj
	dfr0 <- data.frame( "stud.index" = rep(1:RR , CC) , 
				"item.index" = rep(1:CC ,each=RR) , 
				"stud.p" = rep( scores1,CC)  , 
				"item.p" = rep( item.p , each=RR ) , 
				"wgt" = matrix( as.matrix( wgt ) , RR*CC , 1 ) ,
				"p0" = matrix( as.matrix(freq.correct ) , RR*CC , 1 )
						)
#	dfr0$item.plogis <- qlogis( dfr0$item.p)
#	dfr0$stud.plogis <- qlogis( dfr0$stud.p)
	if (progress){ 
		cat("\n*******Logistic Model***********\n") 
#		cat("*** " , paste(Sys.time())  , "*******\n")		
			}	
	theta <- stats::qlogis( scores1 )
	b <- stats::qlogis( item.p )
	#***********************************
	# begin algorithm
    numdiff.parm <- .001
	max.increment <- 1
	deviation <- 1000
	iter <- 0
	# adjust freq.correct if necessary
	adj <- min( c( min( freq.correct[ freq.correct > 0 ] ) , 
			1-max( freq.correct[ freq.correct < 1 ] )	) ) / 3
	freq.correct1 <- freq.correct
	freq.correct1[ freq.correct == 0 ] <- adj
	freq.correct1[ freq.correct == 1 ] <- 1-adj	
	while( ( iter < maxit) & ( deviation > conv ) ){   
		b0 <- b
		theta0 <- theta
		# update theta
		res <- .update.theta.logistic( theta , b , freq.correct=freq.correct1 , wgt , 
					numdiff.parm , max.increment= max.increment)
		theta <- res$theta
		ll <- res$ll
		# update b
		res <- .update.b.logistic( theta , b , freq.correct=freq.correct1 , wgt , 
					numdiff.parm , max.increment= max.increment)
		b <- res$b
		b <- b - mean(b)
	    deviation <- max( abs( c( theta-theta0 , b - b0 )) )
		if (progress){ 
				cat( "Iteration" , iter , "- Deviation =" ,  round( deviation , 6 ) , "\n")
				flush.console()
					}	
        iter <- iter + 1					
			}  # end algorithm
	#*************************
	dfr0$irtfitted <- matrix( plogis( outer(theta , b , "+") ) , RR*CC , 1 )[,1]
	Y <- 0*freq.correct
	Y <- as.matrix(Y)
	Y[ as.matrix(dfr0[ , c("stud.index" , "item.index" ) ] ) ] <- dfr0$irtfitted
	G <- Y
	# fitted frequency table
	dfr0$stud.p <- rep( theta ,CC)
	dfr0$item.p <- rep( b , each=RR )
	colnames(dfr0)[ colnames(dfr0) == "stud.p" ] <- "person.sc"
	colnames(dfr0)[ colnames(dfr0) == "item.p" ] <- "item.sc"	
	colnames(dfr0)[ colnames(dfr0) == "p0" ] <- "freq"		
	colnames(dfr0)[ colnames(dfr0) == "irtfitted" ] <- "freq.fitted"
    dfr0 <- dfr0[ order( dfr0$stud.index * 1000 + dfr0$item.index ) , ]							
	# deviation criterion
	wgt1 <- ( wgt / colSums( wgt ) ) / ncol(wgt)
	fit <- sqrt( sum( ( freq.correct - G  )^2 * wgt1  ) )
	# calculate likelihood
    ll <- list( 
		"ll.ind" = .calc.ll.isop( freq.correct , wgt , irtfitted =freq.correct )	
					)
	ll$ll.logistic <- .calc.ll.isop( freq.correct , wgt , irtfitted =G )
	NW <- mean( colSums(wgt) )
	ll$llcase.ind <- ll$ll.ind/NW	
	ll$llcase.logistic <- ll$ll.logistic/NW
#	if (progress){ 
#		cat("*** " , paste(Sys.time())  , "*******\n")			
#			}	
	# output
	res <- list( "fX" = G , "ResX" = freq.correct - G , 
				"fit" = fit  , "item.sc"= b  , 
				"person.sc" = theta , "ll"=ll  , "freq.fitted"=dfr0)
    return(res)
	}
####################################################################