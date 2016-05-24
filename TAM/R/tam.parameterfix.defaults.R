
###############################################################
# generate design matrices for parameter fixings

generate.xsi.fixed.estimated <- function( xsi , A ){
	L <- length(xsi)
	xsi.fixed.estimated <- 	cbind( 1:L , xsi )
	rownames(xsi.fixed.estimated) <- dimnames(A)[[3]]
	return(xsi.fixed.estimated)
			}
			
###################################################
generate.B.fixed.estimated <- function(B){			
	dimB <- dim(B)
	I <- dimB[1]
	K <- dimB[2]
	D <- dimB[3]
	B.fixed.estimated <- matrix( 0 , nrow=I*K*D , ncol= 4 )
	vv <- 1
	for (ii in 1:I){
	for (kk in 1:K){
	for (dd in 1:D){
		B.fixed.estimated[vv , 1:4] <- c( ii , kk , dd , B[ ii , kk , dd ] )
		vv <- vv + 1 
				}
			}
		}
	return( B.fixed.estimated  )
			}