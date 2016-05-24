
##################################################################
# application of scale for a list of multiply imputed datasets,
# single datasets or nested multiply imputed datasets
scale_datlist <- function( datlist , orig_var , trafo_var , weights = NULL ,
                    M=0, SD=1 , digits = NULL){
	is_dfr <- FALSE
	is_iL <- FALSE
    is_NIL <- FALSE
	is_ndl <- FALSE
	
	if ( class(datlist) == "nestedImputationList"){
			datlist <- nested.datlist_create( datlist )
			is_iL <- TRUE
							}
	
	if ( class(datlist) == "nested.datlist" ){
			datlist <- nested.datlist2datlist(datlist)
			Nimp <- attr(datlist , "Nimp" )
			is_ndl <- TRUE		
							}

	
	if ( class(datlist) == "imputationList"){
	        datlist0 <- datlist
			datlist <- datlist_create( datlist )
			is_iL <- TRUE
							}
							
	
	#**** processing if datlist is a data frame
	if ( ! ( class(datlist) %in% c("datlist") ) ){
	      is_dfr <- TRUE 
	      datlist0 <- datlist 
          datlist <- list( 1 )
		  datlist[[1]] <- datlist0
		  class(datlist) <- "datlist"
					}
					
    #*** create weights if needed					
    PP <- length(datlist)
    if ( is.null(weights) ){
	    N <- nrow(datlist[[1]])
        weights <- rep(1,N)
                            }
	weights0 <- weights
	#---- compute means and standard deviations
    res <- lapply( datlist , FUN = function(dd){
        # dd <- datlist[[1]]
		if ( is.character(weights0) ){
			weights <- dd[ , weights0 ]
					}
        m1 <- stats::weighted.mean( dd[,orig_var] , w = weights )
        sd1 <- sqrt( Hmisc::wtd.var( dd[,orig_var] , weights ) )
        c(m1,sd1)
            } )
    #---- compute averaged mean and SD
	res <- matrix( unlist(res) ,ncol=2 , byrow=TRUE ) 
    a1 <- colMeans(res)
	#---- create derived variable
    for (pp in 1:PP){
        dd <- datlist[[pp]]
        dd[,trafo_var] <- M + SD * ( dd[,orig_var] - a1[1] ) / a1[2]
        if ( ! is.null(digits) ){
                dd[,trafo_var] <- round( dd[,trafo_var] , digits )
                                }
        datlist[[pp]] <- dd
            }
	#---- output
	if ( is_dfr ){
		datlist <- datlist[[1]]
				 }			
	if ( is_iL ){
        datlist0$imputations <- datlist
		datlist <- datlist0
				}				 
	if ( is_ndl ){
		datlist <- datlist2nested.datlist(datlist=datlist, Nimp=Nimp)	
				}
	if ( is_NIL ){
		datlist <- NestedImputationList(datlist)
				}								
    return(datlist)
                }
######################################################################