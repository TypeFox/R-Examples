##############################################
# modify q-matrix
mcdina.modify.qmatrix <- function( q.matrix , skillclasses){	
	# create new q.matrix
	K <- ncol(q.matrix) - 2
	maxattr <- apply( q.matrix[,-c(1:2) ] , 2 , max )
	qmatrix_mod <- NULL
	q.matrix1 <- q.matrix[,1:2]	
	K1 <- max(maxattr)	
	res <- list( "q.matrix" = q.matrix , "q.matrix0" = NULL , 
			"maxmaxattr" = K1 , "skillclasses" = skillclasses ,
			"skillclasses0" = skillclasses , "qmatrix_mod" = NULL )
	if (K1 > 1 ){	
		m1 <- matrix( 0:K1 , nrow=K1+1 , ncol=K )
		skillclasses <- as.matrix( expand.grid( as.data.frame( m1) ) )
		colnames(skillclasses) <- colnames(q.matrix)[ -c(1:2) ]
		# create modified q-matrix
		for (kk in 1:K){ # kk <- 1
			qmatrix_mod.kk <- data.frame( "attr_index" = kk  , 
				"maxattr" = maxattr[kk] )
			skillclasses <- skillclasses[ skillclasses[,kk] <= maxattr[kk] , ]
			for (zz in 1:(maxattr[kk] ) ){ # 	zz <- 1
				name <- paste0( colnames(q.matrix)[kk+2] , ".L" , zz )
				q.matrix1[ ,  name ] <- 1 * ( q.matrix[ , kk + 2] >= zz )
						}
			qmatrix_mod <- rbind( qmatrix_mod , qmatrix_mod.kk )		
					}
		qmatrix_mod$start <- c(1,cumsum( qmatrix_mod$maxattr)[ - K ] + 1  )
		qmatrix_mod$end <- cumsum( qmatrix_mod$maxattr)
		skillclasses0 <- skillclasses	
		rownames(skillclasses0) <- .matrixstring( skillclasses0 , "P" )
		skillclasses <- as.data.frame(skillclasses)		
		# create modified skillclasses
		for (kk in 1:K){ # kk <- 1
			for (zz in 1:(maxattr[kk] ) ){ # 	zz <- 1
				name <- paste0( colnames(q.matrix)[kk+2] , ".L" , zz )
				skillclasses[ ,  name ] <- 1 * ( skillclasses[ , kk ] >= zz )
						}
					}	
		skillclasses <- skillclasses[ , - c(1:K) ]
		rownames(skillclasses) <- .matrixstring( skillclasses , "P" )	
		res$q.matrix <- q.matrix1
		res$skillclasses <- as.matrix(skillclasses)
		res$skillclasses0 <- skillclasses0
		res$q.matrix0 <- q.matrix
		res$qmatrix_mod <- qmatrix_mod		
		}
	return(res)
	}	
########################################################	
	
######################################################
.mcdina.prepare.qmatrix <- function( dat , q.matrix ){	
	if ( min( dat , na.rm=TRUE ) == 0 ){
		dat <- dat + 1
		I <- ncol(dat)
		if ( nrow(q.matrix) == I ){
			q1 <- data.frame( "item" = 1:I , "categ" = 2 , q.matrix )
			q0 <- data.frame( "item" = 1:I , "categ" = 1 , 0+0*q.matrix )			
			q.matrix <- rbind( q0 , q1 )
			q.matrix <- q.matrix[ order( 100 * q.matrix$item + q.matrix$categ ) , ]
					}
				}
	res <- list("dat"=dat , "q.matrix" = q.matrix )
	return(res)
		}


###################################################
# initial estimate of item parameters delta
.mcdina.init.delta <- function( lc , lr ){
    I <- max( lc$item )
	lc$cats <- lc$cats
	#	lc.cats <- lc$cats
    CC <- max(lc$cats)
    delta_ideal <- delta <- array( 0 , dim=c(I , CC , CC ) )
    delta_ideal[ as.matrix( lc[ , c("item" , "cats" , "lr_index" ) ] ) ] <- 1
    eps <- 1E-10
    # define initial delta estimate
    for (ii in 1:I){
	    dii <- delta_ideal[ii,,]
		lcii <- lc[ lc$item == ii , ]
		Ncii <- nrow(lcii)
        for (cc in 1:CC){
            dii.cc <- dii[,cc]
            delta[ii,,cc] <- dii.cc * ( 0.8 / sum( dii.cc + eps) ) +
                                (1-dii.cc) * ( .2 / sum( ( 1-dii.cc) + eps ) )
			if (Ncii < CC ){ 
				delta[ii,seq(Ncii+1,CC),cc] <- 0
				delta[ii,,cc] <- delta[ii,,cc] / sum( delta[ii,,cc] )
							}
            if ( sum( dii.cc ) == 0 ){ delta[ii,,cc] <- 0 }
						
                        }
                    }
    res <- list( "delta" = delta , "delta_ideal" = delta_ideal )
    return(res)
        }
############################################################

##########################################
# preparation function for whole test
.mcdina.prep.test.latent.response <- 
function( q.matrix , K , TP , skillclasses , classes ){
	I <- length( unique(q.matrix[,1]))
	lr <- NULL
	lc <- NULL
	itemstat <- NULL
    for (ii in 1:I){ 
	   res <- .mcdina.prep.item.latent.response( ii , q.matrix , 
				K , TP , skillclasses , classes )
		lr <- rbind( lr , res$lr )
		lc <- rbind( lc , res$lc )
		itemstat <- rbind( itemstat , res$itemstat )
			}
	res <- list("lr"=lr , "lc"=lc , "itemstat"=itemstat)			
	return(res)
		}
###############################################		
		
##############################################
# compute preparation table for one item
.mcdina.prep.item.latent.response <- 
function( ii , q.matrix , K , TP , skillclasses , classes ){
	q.ii <- q.matrix[ q.matrix[,1] == ii , ]	
#	classes <- rownames(skillclasses)
	# categories
	cats.ii <- q.ii[,2]
	CC <- length(cats.ii)
	# calculate relevant attributes
#	qsum <- rowSums( q.ii[ , 1:K + 2  ] )
	qsum <- rowSums( q.ii[ , 1:K + 2  ]  )
	index.max <- which( qsum == max(qsum) )
	# necessary attributes for item ii
	attr.ii <- which( q.ii[ index.max[1] , 1:K + 2] > 0 )
	q.ii.red <- q.ii[ , attr.ii + 2 , drop=FALSE]
	# calculate matrix with skill classes
	sk.ii1 <- sk.ii2 <- matrix( 0 , nrow=TP , ncol=CC)
	colnames(sk.ii1) <- colnames(sk.ii2) <- paste0("Cat" , cats.ii )
	rownames(sk.ii1) <- rownames(sk.ii2) <- rownames(skillclasses)
	for (cc in 1:CC){
		sk.ii2[ , cc] <- 1 * ( rowSums( skillclasses[ , attr.ii , drop=FALSE] != q.ii.red[rep(cc,TP) ,] ) == 0 )
		tmp1 <- skillclasses[ , attr.ii , drop=FALSE] %*% t( q.ii.red[cc,]  )
		sk.ii1[ , cc] <- 1 * ( tmp1 >=  sum( q.ii.red[cc ,] ) ) 
		sk.ii1[ , cc] <-  tmp1*sk.ii1[ , cc]
				}				
	sk.ii1 <- 1 * ( sk.ii1 > 0 )
	v1.ii <- which( rowSums( sk.ii1 ) == 0 )
	i5 <- which( rowSums( q.ii.red ) == 0 )
	sk.ii1[ v1.ii , i5 ] <- 1
	ind.ii <- which( rowSums( sk.ii2 ) == 0 )
	sk.ii2[ind.ii , ] <- sk.ii1[ ind.ii , ]
	
	# define latent response groups
	lg <- "LR"
	for (cc in 1:CC){
		lg <- paste0( lg , ifelse( sk.ii2[,cc]==1 , cats.ii[cc] , "") )
				}
#	sk.ii3 <- cbind( sk.ii2 , lg )
	groups <- sort( unique(lg) )
	lr <- data.frame("item" = ii , "skillclass" = classes , 
		"skillclass_index" = 1:TP , "lr" = lg )
	lr$lr_index <- match( lr$lr , groups )
	# unique latent groups
	lg1 <- sapply( cats.ii , FUN = function(cc){ grep( cc , groups) } )
	lc <- data.frame("item"=ii , "cats"=cats.ii ,
				"lr"= groups[ lg1 ] )
			
	lc$max.cat <- 0
	lc$max.cat[ index.max ] <- 1
	lc$lr_index <- match( lc$lr , groups )
	lc$Q <- .matrixstring( q.ii[ , 1:K + 2  ] , "Q" )
	lc$lr_level <- rowSums( q.ii[ , 1:K + 2  ])
	lc <- lc[ order( paste( lc$lr_level , lc$cats) ) , ]
	lc$lr_level <- paste0( lc$lr_level , 
		LETTERS[ match( lc$lr , unique(lc$lr) ) ] )
	lc <- lc[ order( paste( lc$cats) ) , ]	
	# item statistics
	itemstat <- data.frame("item" = ii , "N.cat" = CC ,
			"N.lr" = length(groups) )
	itemstat$N.attr <- length(attr.ii)
	res <- list("lr"=lr , "lc"=lc , "itemstat"=itemstat)			
	return(res)
		}
############################################################
# calculates a string pattern consisting of matrix entries
# matr <- skillclasses
# string <- "Q"
.matrixstring <- function( matr , string ){
	VV <- ncol(matr)
	l1 <- string
	for ( vv in 1:VV){
		l1 <- paste0( l1 , matr[,vv] )
				}
	return(l1)
		}
#################################################################