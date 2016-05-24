equivalent.dina <-
function( q.matrix , reparametrization="B" ){
    K <- ncol(q.matrix)
    I <- nrow(q.matrix)
    if ( is.null( colnames(q.matrix)) ){ 
        colnames(q.matrix) <- paste( "S" , 1:K , sep="")
                }
    if ( is.null( rownames(q.matrix)) ){ 
        rownames(q.matrix) <- paste( "I" , 1:I , sep="")    
                }
    
    L <- 2^K 
   
    attr.patt <- matrix( rep( 0, K*L) , ncol=K)
    h1 <- 2
	if (K>=2 ){
		for(ll in 1:(K-1) ){
			lk <- utils::combn( 1:K, ll ) 
			lk
			for ( jj in 1:( ncol(lk) ) ){ 
				attr.patt[ h1, lk[,jj] ] <- 1
				h1 <- h1 + 1
				}
			}
		}
    attr.patt[ L, ] <- rep( 1, K )
	alpha <- attr.patt
	l1 <- apply( alpha , 1 , FUN = function( hh){ paste( hh , collapse="") } )

    # define alpha classes
    alpha <- matrix( NA ,  L , K )
    for (kk in 1:K){ alpha[,kk] <- as.numeric(substring(l1,kk,kk)) }
    alpha <- data.frame(alpha)
    rownames(alpha) <- l1
    colnames(alpha) <- colnames(q.matrix)
    
	#########################################################
	# reparametrization B
	if (reparametrization=="B" ){	
		# equivalent q.matrix: q.matrix*
		qclasses <- apply(q.matrix , 1 , FUN = function(ll){ paste(ll, collapse="" ) } )
		q.matrix.ast <- matrix( 0 , I , L )
		rownames(q.matrix.ast) <- rownames(q.matrix)
		colnames(q.matrix.ast) <- paste( "S*" , rownames(alpha) , sep="")
		q.matrix.ast[ cbind( 1:I , match( qclasses , paste(rownames(alpha )) ) ) ] <- 1
		q.matrix.ast <- q.matrix.ast[ , - 1 ]
		
		# equivalent alpha: alpha*
		alpha.ast <- matrix( 0 , L , L  )
		rownames(alpha.ast) <- rownames(alpha)
		colnames(alpha.ast) <- rownames(alpha)
		for (tt in 1:L ){ 
			# tt <- 2     # class tt
			alpha.tt <- alpha[ tt , ]    
			q.tt <- alpha^( alpha.tt[ rep(1,L) , ] )
			q.tt <- apply( q.tt , 1 , prod )
			names(q.tt) <- rownames(alpha)
			r.tt <- rowsum( q.tt , names(q.tt) )
			ind.tt <- rownames(r.tt)[ r.tt[,1] > 0  ]
			alpha.ast[ tt , ind.tt ] <- 1
				}
		alpha.ast <- ( t( alpha.ast ) )[,-1]
			}
	######################################################
	# reparametrization A
	if (reparametrization=="A" ){	

		# alpha*
		alpha.ast <- matrix( 0 , L , L )
		rownames(alpha.ast) <- rownames(alpha)
		colnames(alpha.ast) <- rownames(alpha)
		diag(alpha.ast) <- 1
		alpha.ast <- alpha.ast[,-1]
		
		# q.matrix*
	 q.matrix.ast <- matrix(0 , I , L )
	 rownames(q.matrix.ast) <- rownames(q.matrix)
	 colnames(q.matrix.ast) <- paste( "S*" , rownames(alpha) , sep="")
	 for (ii in 1:I){
		#ii <- 2
		q.ii <- q.matrix[ii,]
		a1 <- apply( alpha^( matrix( q.ii , nrow=L , ncol=K , byrow=TRUE )) , 1 , prod )
		q.matrix.ast[ii,] <- a1
				}
		q.matrix.ast <- q.matrix.ast[,-1]
		}
	
	#*****************************
	# collect results
    res <- list( "q.matrix" = q.matrix  , "q.matrix.ast" = q.matrix.ast )
    res$alpha <- alpha
    res$alpha.ast <- alpha.ast
    return(res)
        }
