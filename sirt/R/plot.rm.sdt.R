

##################################################
# plot.rm.sdt
plot.rm.sdt <- function( x , ask=TRUE , ... ){
	rt1 <- x$rater2
	items <- unique( rt1$item )
	maxK <- x$maxK
	for (ii in items){
		# ii <- items[1]
		rtii <- rt1[ rt1$item == ii , ]
		rtii1 <- rtii[ , grep( ".trans" , colnames(rtii) ) ]
		# xlim <- range(rtii1)
		# xlim <- c(-1, ncol(rtii1) + 1 )
		xlim <- c(-1 , maxK[ii] + 1 )
		# K <- ncol(rtii1)
		K <- maxK[ii]
		for (kk in 1:K){
			# kk <- 1
			m1 <- as.vector(rtii1[,kk])
			if (kk==1){
				graphics::dotchart( m1 , xlim=xlim , labels= rtii$rater , 
					main = paste0("Rater Parameters | Item " , ii ) , ...)
						}
			RR <- nrow(rtii)
			graphics::lines( rep( kk - .5 , RR ) , 1:RR , lty=2 , col= "darkred" )
			graphics::points( m1 , 1:RR , pch=16 )
			graphics::lines( m1 , 1:RR , col="darkred" )    
						}
		graphics::text( K + 1 - 0.2 , 1:RR , format( rtii$d , digits=2, nsmall=2 )	)
		graphics::text( K + 1 - 0.2 , RR + .7  , expression(d) )
		graphics::par( ask=ask )
			}
		}
###########################################################		

