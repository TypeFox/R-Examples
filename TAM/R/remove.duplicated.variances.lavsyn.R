

			
#*****************************************************			
# remove duplicated variances
remove.duplicated.variances.lavsyn <- function( res0 , items ) {
	res0 <- gsub( " " , "" , res0 )
	res1 <- strsplit( res0 , split="\n")[[1]]
	res1 <- data.frame( "syn" = res1 , "sel" = 0 )
	res1$variance <- 0
	ind <- grep( "~~" , res1$syn )
	res1$variance.index <- 0
	l0 <- strsplit( paste(res1[ ind, "syn" ]) , split="~~")
	l0 <- unlist( lapply( l0 , FUN = function(ll){ ll[1] } ) )
	ind <- ind[ l0 %in% items ]

	if ( length(ind) > 0){
		res1$variance.index[ind] <- 1
					}
	
	if ( length(ind) > 0 ){
		res1[ ind,"variance"] <- 1
		res1$variance.obs <- ""
		l1 <- res1[ ind, "syn" ]
		l1 <- strsplit( paste(l1) , split="~~")
		l2 <- lapply( l1 , FUN = function(ll){ ll[1] } )
		res1$variance.obs[ind] <- unlist(l2)
		l3 <- duplicated( res1[ind , "variance.obs"] )
		if ( sum(l3) > 0 ){
			res1 <- res1[ - ind[ l3 ] , ]
							}
					}
	# recreate  lavaan syntax
	lav2 <- paste( res1$syn , collapse="\n") 			
	return(lav2)
			}

