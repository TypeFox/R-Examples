

###########################################
# Rcpp version of R's table function
fasttable <- function( vec  , sort.names=FALSE ){
    datavec <- matrix( vec , ncol=1 )
    # res <- bifie_fasttable( datavec )	
	if ( storage.mode(vec) == "character" ){
				characters <- TRUE	
					} else {
				characters <- FALSE
						}
	if ( ! characters ){
		res <- .Call("bifie_fasttable" , datavec , PACKAGE="BIFIEsurvey" )
		res1 <- res$tableM[ 1:res$N_unique , , drop=FALSE]
		tvec <- res1[,2]
		names(tvec) <- res1[,1]
				}
	if ( characters ){ 			
		# t1 <- bifie_table1_character( vec )
		t1 <- .Call("bifie_table1_character" , vec , PACKAGE="BIFIEsurvey")
		res <- t1$tableM
		names(res) <- t1$table_names
		if ( sort.names ){ 
			tvec <- res[ sort( names(res) ) ]
					} else { tvec <- res }
					}								
    return(tvec)
        }
###########################################		

bifietable <- fasttable