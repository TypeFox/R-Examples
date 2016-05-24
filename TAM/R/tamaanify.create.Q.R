


#############################################################
# create Q matrix
tamaanify.create.Q <- function( res ){
	lavpartable <- res$lavpartable
	resp <- res$resp
	ind <- which( lavpartable$op == "=~" )
	#*** set some values in ustart
	lv1 <- lavpartable[ ind , ]
	lv1[ is.na(lv1$ustart) , "ustart" ] <- 1
	lavpartable[ ind , ] <- lv1 
	
	facs <- unique( paste( lavpartable$lhs[ind] )) 
	NF <- length(facs)
	items <- colnames(resp)
	I <- length(items)
	Q <- matrix( 0 , nrow=I , ncol=NF)
	rownames(Q) <- items
	colnames(Q) <- facs
	for (ff in 1:NF){
		# ff <- 1
		ind.ff <- intersect(   which( paste(lavpartable$lhs) == facs[ff] ) , ind )
		lff <- lavpartable[ ind.ff , ]	
		Q[ rownames(Q) %in% paste(lff$rhs) , ff ] <- lff$ustart
					}
	res$Q <- Q					
	res$lavpartable <- lavpartable
	# res <- list( "Q" = Q , "lavpartable" = lavpartable )
	return(res)
		}
#############################################################