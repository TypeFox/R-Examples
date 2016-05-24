


##################################################
# create table of item parameters
.TAM.itempartable <- function( resp , maxK , AXsi , B , ndim ,
			resp.ind , rprobs,n.ik,pi.k){
	
	if ( is.null(dimnames(B)[[1]] ) ){ 
		dimnames(B)[[1]] <- colnames(resp)
					}
		
	item1 <- data.frame( "item" = dimnames(B)[[1]] )
	item1$N <- colSums(resp.ind )
	item1$M <- colSums( resp.ind * resp , na.rm=TRUE) / colSums( resp.ind )
	
	maxKi <- rowSums( 1 - is.na( AXsi ) ) - 1
	I <- nrow(item1)

	item1$xsi.item <- - AXsi[ cbind(1:I , maxKi+1) ] / maxKi
#	item1$xsi.item <- AXsi[ cbind(1:I , maxKi+1) ] 

	
	#****
	# Item fit
	# probs ... [ classes , items , categories ]
	probs <- aperm( rprobs , perm=c(3,1,2))
	pi.k <- matrix( pi.k , ncol=1 )
#	res <- .tam.itemfit.rmsea( n.ik , pi.k , probs )
	#####
	# Exploratory analyses show that item fit rmsea
	# does not seem to be sensitive
#	item1$rmsea <- res

	b0 <- sum( B[ , 1 , ] , na.rm=TRUE ) 
	# a0 <- sum( A[ , 1 , ] , na.rm=TRUE ) 
	a0 <- 0
	if ( b0 + a0 > 0 ){
		kvec <- 0:(maxK-1)		
				} else {
		kvec <- 1:(maxK-1)
						}
	for (kk in kvec){ # kk <- 1
		item1[ , paste0("AXsi_.Cat" , kk) ] <- - AXsi[,kk+1]
						}
	for (kk in kvec){ # kk <- 1
		for (dd in 1:ndim){
			item1[ , paste0("B.Cat" , kk,".Dim",dd) ] <- B[,kk+1,dd]
							}
					}						
    item1 <- item1[ item1$N > 0 , ]	
#	item1 <- item1[ order( paste( item1$item)) , ]		
	rownames(item1) <- NULL
	
	#*** substitute -99 by missing
	item1[ item1 == -99 ] <- NA
	
	return(item1)
		}