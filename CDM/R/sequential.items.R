#################################################
# converts a dataset with polytomous item responses
# into a dataset with sequential dichotomous items
sequential.items <- function( data ){
	N <- nrow(data)
	I <- ncol(data)
	maxK <- apply( data , 2 , max , na.rm=TRUE )
	dat.exp  <- matrix( NA , nrow=N , ncol=sum(maxK) )
	vv <- 1
	for (ii in 1:I){
		# ii <- 4     # item
		dat.ii <- data[,ii]    
		kk <- 1
		dat.exp[ , vv ] <- 1 * ( dat.ii >= kk )
		vv <- vv + 1     # column index
		if (maxK[ii]>1){
			for (kk in 2:maxK[ii]){
				dat.exp[ , vv ] <- 1 * ( dat.ii >= kk )
				dat.exp[ dat.ii < kk - 1 , vv ] <- NA
				vv <- vv + 1            
							}
				}        
			}
	#****
	# variable names
	varnames <- sapply( 1:I , FUN = function(ii){
			if ( maxK[ii] == 1 ){  v1 <- colnames(data)[ii] } else
				{ v1 <- paste0( colnames(data)[ii] , "_Cat" , 1:maxK[ii] ) }
			v1
			  } )
	colnames(dat.exp) <- unlist( varnames  )
	dat.exp <- as.data.frame( dat.exp)
	# item information table
	iteminfo <- data.frame("item" = rep( colnames(data) , maxK ) )
	iteminfo$itemindex <- match( iteminfo$item , colnames(data) )
	iteminfo$category <- unlist( sapply( 1:I , FUN = function(ii){ 1:(maxK[ii])  } , simplify=FALSE)   )
	iteminfo$pseudoitem <- colnames(dat.exp)
	res <- list( "dat.expand" = dat.exp , "iteminfo" = iteminfo ,
			"maxK"= maxK)
	return(res)
		
		}
##############################################################