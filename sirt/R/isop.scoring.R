#################################################
# scoring students and items according to the
# Scheiblechner's ISOP model
isop.scoring <- function( dat , score.itemcat=NULL){
	#***
	# compute modified percentile score
	maxK <- apply( dat , 2 , max , na.rm=TRUE )
	K <- max(maxK)
	I <- ncol(dat)
	# define response matrix
	dat.resp <- 1-is.na( dat )
	# p values for every category
	p.itemcat <- matrix(NA , nrow=I , ncol=K+1 )
	rownames(p.itemcat) <- colnames(dat)
	colnames(p.itemcat) <- paste0( "p.Cat" , 0:K )
	for (kk in 0:K){
		p.itemcat[,kk+1] <- colMeans( dat==kk , na.rm=TRUE )
				}				
	#**************************************				
	if ( is.null( score.itemcat )){
		score.itemcat <- 0*p.itemcat					
		colnames(score.itemcat) <- paste0( "score.Cat" , 0:K )	
		# score itemcat 
		for (kk in 0:K){
			if (kk>0){
				v1 <- apply( p.itemcat[ , 1+seq( 0 , kk -1 , 1 ) , drop=FALSE] ,
							1 , sum , na.rm=TRUE )
					} else { v1 <- 0 }
			if (kk<K){
				v2 <- apply( p.itemcat[ , 1+seq( kk+1 , K , 1 ) , drop=FALSE] ,
							1 , sum , na.rm=TRUE )
					} else { v2 <- 0 }	
			score.itemcat[,kk+1] <- v1-v2
					}
				}  # end score itemcat
	  #*************************								
	  dat <- dat[ , rownames(score.itemcat) ]
	  # scoring table
	  score.dat <- 0*dat
	  for (kk in 0:K){
		score.dat <- score.dat + matrix( score.itemcat[,kk+1] , nrow=nrow(dat) , 
						ncol=I , byrow=TRUE ) * ( dat == kk ) 
					}
	  # score persons
	  person <- data.frame( "pid" = 1:nrow(dat) ,
	                "score" = rowSums( dat , na.rm=TRUE ) ,
					"max" = rowSums( dat.resp ) , 
					"M" = rowMeans( dat , na.rm=TRUE ) ,
					"mpsc" = rowMeans( score.dat , na.rm=TRUE )
							)
	  # compute item p score
	  item.pscore <- rep(NA,I)
	  for (ii in 1:I){
#		  ii <- 1
		  v1 <- colMeans( dat[,ii ] > dat[,-ii] )
		  v2 <- colMeans( dat[,ii ] < dat[,-ii] )
		  item.pscore[ii] <- mean( v1 - v2 )
						}
	  # score items
	  item <- data.frame( "item" = colnames(dat) ,
			   "M" = colMeans( dat , na.rm=TRUE ) )
	  item$pscore <- item.pscore
	  item <- cbind( item , p.itemcat ,  score.itemcat )
	  # calculate distribution function
	  distr.fct <- t(apply( p.itemcat , 1 , FUN = function(ll){
					cumsum( ll ) } ) )									
#	  colnames(item)[-c(1:2)] <- colnames( score.itemcat )[-1]
	  #***
	  # collect results
	  res <- list( "person"=person , "item"=item , 
			"p.itemcat"=p.itemcat , "score.itemcat"=score.itemcat ,
			"distr.fct"=distr.fct)
	  return(res)
	  }
##############################################################
# example scoring for item with 3 categories
# probabilities p0 , p1 , p2
# score category 0: s0 = - ( p1 + p2 )
# score category 1: s1 = p0 - p2
# score category 2: s2 = p0 + p1
# ***
# s2 + s0 = p0 + p1 - p1 - p2 = p0 - p2 = s1
# ***
# s2 - s1 = p0 + p1 - p0 + p2 = p1 + p2
# s1 - s0 = p0 - p2 + p1 + p2 = p0 + p1
