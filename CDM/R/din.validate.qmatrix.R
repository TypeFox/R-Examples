######################################################
# Q-matrix validation
din.validate.qmatrix <- function( object , digits = 3 , print=TRUE){
	#*****
	# extract original Q-matrix
	q.matrix <- object$q.matrix
	rule <- object$rule
	# extract estimated parameters
	guess <- object$guess[,1]
    slip <- object$slip[,1]
	# calculate originally estimated IDI
	IDI <- 1-slip-guess
	#****
	# define all possible Q-matrix vectors
	nodes <- c(0,1)
	K <- ncol(q.matrix)
	L <- 2^K
	q.matrix.poss <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, K) , ncol = K ) ) ) )	
	colnames(q.matrix.poss) <- colnames(q.matrix)
	q.matrix.poss <- q.matrix.poss[ ! ( rowMeans( q.matrix.poss ) %in% c(0) ) , ]
	QQM <- nrow(q.matrix.poss) 	# number of possible Q-matrix vectors
	#****
	# extract data and attributes
	data <- object$data	
	I <- ncol(data)		# number of items
	I.lj <- object$I.lj
	R.lj <- object$R.lj
	attr.patt <- object$attribute.patt.splitted
	# calculate modification item parameters
	coef.modified <- matrix( 0 , nrow=QQM*I , ncol=4 )
	colnames(coef.modified) <- c("item" , "qmatrix.row" , "guess" , "slip" )
	coef.modified <- as.data.frame( coef.modified )
	coef.modified$item <- rep( 1:I , QQM )
	coef.modified$qmatrix.row <- rep( 1:QQM , each=I )	
	for (qqm in 1:QQM){	
#		qqm <- 1
		qt.matrix <- q.matrix.poss[ rep(qqm,I) , ]
#		compt <- rowSums( qt.matrix )	
		compt <- ( rowSums(qt.matrix)  )*( rule=="DINA")   + 1* ( rule == "DINO" )  		
		# calculate expected counts
		I.j0 <- rowSums( ( t(( attr.patt %*% t(qt.matrix) ) )  <  outer(  compt, rep(1,L) ) ) * I.lj ) 
		I.j1 <- rowSums( ( t(( attr.patt %*% t(qt.matrix) ) )  >= outer(  compt, rep(1,L) ) ) * I.lj )             
		R.j0 <- rowSums( ( t(( attr.patt %*% t(qt.matrix) ) )  <  outer(  compt, rep(1,L) ) ) * R.lj )
		R.j1 <- rowSums( ( t(( attr.patt %*% t(qt.matrix) ) )  >= outer(  compt, rep(1,L) ))  * R.lj )	
		# calculate guessing and slipping parameters
#		I.j0[ I.j0 == 0] <- 0.05
#		I.j1[ I.j1 == 0] <- 0.05           
		guess <- R.j0 / I.j0
		slip <- ( I.j1 - R.j1 ) / I.j1
		coef.modified[ I*( qqm - 1 ) + 1:I , "guess" ] <- guess
		coef.modified[ I*( qqm - 1 ) + 1:I , "slip" ] <- slip
						}
	coef.modified <- coef.modified[ order( coef.modified$item ) , ]
	coef.modified$IDI <- 1 - coef.modified$slip - coef.modified$guess
	# look for original rows
	coef.modified$qmatrix.orig <- 1 * ( rowMeans( q.matrix.poss[ coef.modified$qmatrix.row , ] 
				==	q.matrix[ coef.modified$item , ] ) == 1 )
	coef.modified$IDI.orig <- IDI[ coef.modified$item ]
	coef.modified$delta.IDI <- coef.modified$IDI - coef.modified$IDI.orig
	# restructure matrix coef.modified
	coef.modified <- data.frame( "item" = colnames(data)[ coef.modified$item ] , 
				"itemindex" = coef.modified$item ,
				q.matrix.poss[ coef.modified$qmatrix.row , ] ,
				coef.modified[ , - c(1:2) ]
						)
    coef.modified[ , -1 ] <- round( coef.modified[,-1] , digits )
	# calculate maximum delta index per item
	a1 <- stats::aggregate( coef.modified$IDI , list( coef.modified$itemindex ) , max )
	coef.modified$max.IDI <- a1[ coef.modified$itemindex , 2]
	coef.modified <- coef.modified[ order( coef.modified$itemindex - coef.modified$IDI ) , ]	
	# print output
	coef.modified2 <- coef.modified
	coef.modified2 <- coef.modified2[ coef.modified2$IDI > coef.modified$IDI.orig , ]
	nochange <- nrow(coef.modified2) == 0
	# calculate proposed Q-matrix
	q.matrix.prop <- q.matrix
    if ( ! nochange ){ 
		items <- unique( coef.modified2$itemindex )
		for (ii in items ){
		    c2 <- ( coef.modified2[ coef.modified2$itemindex == ii , ] )[1,]
			q.matrix.prop[ ii , ] <- as.vector( t(c2[ 1, seq(3 , 3+K -1) ]) )
						}
				}	
	if ( print ){
		if ( ! nochange ){ 		
				print( coef.modified2 ) 
				cat("\nProposed Q-matrix:\n\n")
				print(q.matrix.prop)
					}
		if ( nochange ){ cat("No Q-matrix entries should be changed.\n") }
				}	
	res <- list( "coef.modified" = coef.modified ,
				"coef.modified.short" = coef.modified2 ,
				"q.matrix.prop" = q.matrix.prop	)
	invisible(res)
	}
######################################################