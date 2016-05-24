



##########################################################
# variance constraints
tamaanify.variance.fixed <- function( res ){
	lavpartable  <- res$lavpartable
	Q <- res$Q
	variance.fixed <- NULL
	facs <- colnames(Q)
	ind1 <- which( paste(lavpartable$lhs) %in% facs )
	ind2 <- which( paste(lavpartable$rhs) %in% facs )
	ind3 <- which( paste(lavpartable$op) == "~~" )
	ind4 <- which( lavpartable$free == 0 )
	ind <- intersect( intersect( intersect(ind1,ind2) , ind3 ) , ind4 )
	lv1 <- lavpartable[ind,]
	LI <- length(ind)
	if (LI>0){
		for ( ii in 1:LI ){
			#  ii <- 3
			v1 <- cbind( match( lv1[ii,"lhs"] , facs ) , match( lv1[ii,"rhs"] , facs ) ,
								lv1[ii,"ustart"] )
			variance.fixed <- rbind( variance.fixed , v1 )
								}
						}
   res$variance.fixed <- variance.fixed
   return(res)
				}
##########################################################
