getOddsRatio <-
function(x, za=TRUE){
	xx <- x
	if( any( x == 0 ) ){xx <- x + 0.5} # agresti2002, p.71
	xx <- log(xx)
	or <- NA

	chk1 <- chk2 <- FALSE
	if( length(dim(x)) == 2 ){
		or <- exp( ( xx[1,1] + xx[2,2]) - ( xx[1,2] + xx[2,1] ) )
		if(za == FALSE){
			chk1 <- any(c(x[1,1], x[2,2] )==0)
			chk2 <- any(c(x[2,1], x[1,2] )==0)
		}
	}
	if( length(dim(x)) == 3 ){ #Bartlett's ratio of odds ratios
		or <- exp( ( xx[1,1,1] + xx[2,2,1] + xx[1,2,2] + xx[2,1,2] ) - ( xx[1,2,1] + xx[2,2,2] + xx[1,1,2] + xx[2,1,1] ) )
		if(za == FALSE){
			chk1 <- any(c(x[1,1,1], x[2,2,1], x[1,2,2], x[2,1,2] )==0)
			chk2 <- any(c(x[1,2,1], x[2,2,2], x[1,1,2], x[2,1,1] )==0)
		}
	}

	if(za==FALSE){if( chk1 & chk2 ){ or <- NaN }else{ or <- ifelse( chk1, 0, Inf ) }}

return(or)
}
