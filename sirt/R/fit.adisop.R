###################################################################
# Fit ADISOP model
fit.adisop <- function( freq.correct , wgt , conv = .0001 , 
		maxit=100 , epsilon=.01 , progress=TRUE , calc.ll=TRUE ){
	#########################
	eps2 <- epsilon
	M1 <- as.matrix( freq.correct )
	wgt <- as.matrix(wgt)
	freq.correct1 <-  ( freq.correct + eps2/2 ) / ( 1 + eps2 )
	wgt2 <- wgt / sqrt( ( freq.correct1 * ( 1 - freq.correct1 ) ) )
	RR <- nrow(freq.correct)
	CC <- ncol(freq.correct)
	# auxiliary variables
	scores <- 0:(RR-1)
	item.psx <- colSums( freq.correct * wgt ) / colSums( wgt )
	# initialization
	dfr0 <- data.frame( "stud.index" = rep(1:RR , CC) , 
				"item.index" = rep(1:CC ,each=RR) , "stud.p" = rep( scores,CC)  , 
				"item.p" = rep( item.psx , each=RR ) , 
				"wgt2" = matrix( as.matrix( wgt2 ) , RR*CC , 1 ) ,
				"p0" = matrix( as.matrix(freq.correct ) , RR*CC , 1 )
						)
	dfr0$studPitem.p <- dfr0$stud.p + dfr0$item.p                    
	dfr0 <- dfr0[ order( dfr0$p0) , ]
#	dfr0$psort.index <- seq(1 , nrow(dfr0)) #  / nrow(dfr0) 
	dfr0$psort.index <- base::rank( dfr0$p0 )
	# define matrix I
	I.matrix <- 0*M1
	I.matrix[ as.matrix(dfr0[ , c("stud.index" , "item.index" ) ] ) ] <- 
				dfr0$psort.index
	IS <- sum(I.matrix)
	Y <- Y0 <- I.matrix
	iter <- 0
	deviation <- 1000
	if (progress){ cat("\n*******ADISOP Model*********\n") }
	#########################
	# begin algorithm
	while( ( iter < maxit) & ( deviation > conv ) ){   
		#****
		# algorithm
		Y0 <- Y
		# row and column averages
#	    X1 <- rowSums(Y0*wgt2) / rowSums( wgt2 )
		X1 <- rowMeans(Y0)
		X1 <- X1 - X1[1] 
#	    X2 <- colSums(Y0*wgt2) / colSums( wgt2 )
		X2 <- colMeans(Y0) 
		X2 <- X2 - X2[1] 
		# calculate sum
		Y <- outer( X1 , X2 , "+" )
		# preparation isotonic regression
		dfr0 <- data.frame( "stud.index" = rep(1:RR , CC) , 
					"item.index" = rep(1:CC ,each=RR) , 
					"stud.X" = rep( X1,CC)  , 
					"item.X" = rep( X2 , each=RR ) , 
					"wgt2" = matrix( as.matrix( wgt2 ) , RR*CC , 1 ) ,
					"Y" = matrix( as.matrix( Y) , RR*CC , 1 ) ,
					"I" = matrix( as.matrix( I.matrix) , RR*CC , 1 ) 
							)    
			# Y = X1 + X2
		dfr0 <- dfr0[ order( dfr0$I) , ]
#		dfr0$yf <- monoreg( x=dfr0$Y , w=dfr0$wgt2 )$yf    
		dfr0$yf <- monoreg.rowwise( matrix( dfr0$Y , nrow=1 ) , matrix(dfr0$wgt2,nrow=1) )[1,]
		Y <- 0*M1
		Y[ as.matrix(dfr0[ , c("stud.index" , "item.index" ) ] ) ] <- dfr0$yf
		Y <- Y / sum(Y) * IS
		# calculate deviation
		deviation <- max( abs( Y - Y0 ))
		iter <- iter + 1
		if (progress){ 
				cat( "Iteration" , iter , "- Deviation =" ,  round( deviation , 6 ) , "\n")
				utils::flush.console()
					}					
				}  # end algorithm 
    #########################				
    #*****
	# calculate link function
	dfr0 <- data.frame( "stud.index" = rep(1:RR , CC) , 
				"item.index" = rep(1:CC ,each=RR) , 
				"stud.p" = rep( X1,CC)  , 
				"item.p" = rep( X2 , each=RR ) , 
				"wgt" = matrix( as.matrix( wgt ) , RR*CC , 1 ),
				"wgt2" = matrix( as.matrix( wgt2 ) , RR*CC , 1 ) ,
				"p0" = matrix( as.matrix(freq.correct ) , RR*CC , 1 )
						)
	dfr0$X1PX2 <- dfr0$stud.p + dfr0$item.p       
	dfr0 <- dfr0[ order( dfr0$X1PX2) , ]           
	dfr0$yf <- monoreg.rowwise( matrix( dfr0$p0 , nrow=1 ) , matrix(dfr0$wgt2,nrow=1) )[1,]	
#	dfr0$yf <- monoreg( x=dfr0$p0 , w=dfr0$wgt2 )$yf
    dfr0 <- dfr0[ order( dfr0$stud.index * 1000 + dfr0$item.index ) , ]							
	
	Y <- 0*M1
	Y[ as.matrix(dfr0[ , c("stud.index" , "item.index" ) ] ) ] <- dfr0$yf
	G <- Y	
	colnames(dfr0)[ colnames(dfr0) == "stud.p" ] <- "person.sc"
	colnames(dfr0)[ colnames(dfr0) == "item.p" ] <- "item.sc"	
	colnames(dfr0)[ colnames(dfr0) == "p0" ] <- "freq"		
	colnames(dfr0)[ colnames(dfr0) == "yf" ] <- "freq.fitted"				
	# deviation criterion
	wgt1 <- ( wgt / colSums( wgt ) ) / ncol(wgt)
	fit <- sqrt( sum( ( M1-G  )^2 * wgt1  ) )
	#****
	# calculate likelihood
	ll <- NULL
	if (calc.ll){
		ll <- list( 
			"ll.ind" = .calc.ll.isop( freq.correct , wgt , irtfitted =freq.correct )	
						)
		ll$ll.adisop <- .calc.ll.isop( freq.correct , wgt , irtfitted =G )
		NW <- mean( colSums(wgt) )
		ll$llcase.ind <- ll$ll.ind/NW	
		ll$llcase.adisop <- ll$ll.adisop/NW
			}
	# output
	res <- list( "fX" = G , "ResX" = M1 - G , "fit" = fit  , "item.sc"= X2 , 
				"person.sc" = X1 , "ll"=ll , "freq.fitted"=dfr0)
	return(res)
	}
############################################################################