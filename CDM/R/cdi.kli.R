
#######################################################
# cognitive diagnostic indices
cdi.kli <- function( object ){
	# object must be of class din or gdina
	if ( ! ( class(object) %in% c("din","gdina") ) ){
		stop("This functions only supports objects of class din or gdina!")
					}
	items <- colnames( object$data )				
	q.matrix <- object$q.matrix
	pjk0 <- pjk <- object$pjk       # [ items , categories , skills ]
	skillclasses <- as.matrix( object$attribute.patt.splitted )
	# rearrange probabilities
	pjk <- aperm( pjk , c(1,3,2 ) )
	dpjk <- dim(pjk)
	pjk <- matrix( pjk , nrow=dpjk[1] , ncol=dpjk[2]*dpjk[3] )

	eps <- 10^(-7) 	# prevent division by zero
	pjk <- ( pjk + eps ) / ( 1 + 2*eps )

	#*****
	# apply Rcpp function for calculation
#	res0 <-  cdm_kli_id( pjk , skillclasses )
	res0 <- .Call( "cdm_kli_id_C", 
					   pjk , skillclasses , 
					   PACKAGE = "CDM")	
	
	#****
	# arrange Kullback Leibler information
	kli <- array( res0$kli , dim=c( res0$TP , res0$TP , res0$I ) )
	kli <- aperm( kli , c(3,1,2) )

	# arrange output
	res <- list( "test_disc" =  sum(res0$glob_item) , "attr_disc" = colSums( res0$attr_item ) , 
				"glob_item_disc" = res0$glob_item , "attr_item_disc" = res0$attr_item ,
				"KLI" = kli , 
				"skillclasses" = res0$skillclasses , "hdist" = res0$hdist , "pjk" = pjk0 ,
				"q.matrix" = q.matrix )					
	# complete summary in a table
    dfr <- cbind( res0$glob_item , res0$attr_item ) 	
	l1 <- c( sum(res0$glob_item) , colSums( res0$attr_item ) )
    dfr <- rbind( l1 , dfr )
	rownames(dfr) <- NULL	
	colnames(dfr) <- c( "cdi_test" , paste0( "cdi_skill" , 1:( ncol(dfr) -1 ) ) )
	dfr <- data.frame( "item" = c("test" , items ) , dfr )
	# names
	names(res$attr_disc) <- colnames(res$attr_item_disc) <- colnames(q.matrix)
	dimnames(res$KLI)[[1]] <- items
	names(res$glob_item_disc) <- rownames(res$attr_item_disc) <- items				
	res$summary <- dfr
	class(res) <- "cdi.kli"
	return(res)
	}
##############################################################

#################################################################################		
# summary S3 method
summary.cdi.kli <- function( object , digits=2, ...){
    obji <- object$summary
    V <- ncol(obji)
    for (vv in 2:V){
        obji[,vv] <- round( obji[,vv] , digits)
                    }
    rownames(obji) <- NULL
    print(obji)
        }
#####################################################################################