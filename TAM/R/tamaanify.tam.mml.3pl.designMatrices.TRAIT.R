
######################################
# TRAIT
tamaanify.tam.mml.3pl.designMatrices.TRAIT <- function( res ){
	anlist <- res$ANALYSIS.list
	items <- colnames(res$resp)
	I <- length(items)
	itemtable <- res$items

	#*****
	# gammaslope parameters
	# res$gammaslope.fixed <- NULL
	res$gammaslope.prior <- NULL
	res$gammaslope.des <- "2PL"
	
	# define appropriate E matrix
	# items x categories x dimensions x parameters
	
	# define also a B.fixed input option
	Q <- res$Q
	Q.fixed <- NULL
	B_fix <- res$B.fixed
	if ( ! is.null(B_fix) ){
		Q.fixed <- NA*Q
		colnames(B_fix) <- c("item_index" , "cat" , "dim" , "value") 	
		B_fix <- B_fix[ B_fix[,"cat"] == 2 , ]
		h1 <-  B_fix[ , c("item_index" , "dim") ]
		Q.fixed[ h1 ] <- B_fix[ , "value"]
						}
	res$Q.fixed <- Q.fixed					
						
	#*****
	# calculate guessing parameters	
	res$est.guess <- NULL
	res$guess <- NULL
	res$guess.prior <- NULL
	lavpartable <- res$lavpartable
	ind <- which( lavpartable$op == "?=" )
	lav1 <- lavpartable[ ind , ]
		
	if ( nrow(lav1) > 0){
		# include labels
		labs <- paste0( lav1$lhs , "_guess" )
		lav1$label <- ifelse( paste(lav1$label) != "" , 
				paste(lav1$label) , labs )
		lavpartable[ind,] <- lav1
		res$lavpartable <- lavpartable
		I <- length(items)
		est.guess <- rep(0,I)
		names(est.guess) <- items
		lav1$item.index <- match( paste(lav1$lhs) , items )
		lav1$ustart <- ifelse( is.na( lav1$ustart) , .2 , lav1$ustart )						
		guess.labels <- unique(paste(lav1$label))
					
		est.guess[ lav1$item.index] <-
				match( paste(lav1$label) , guess.labels )
		guess <- 0 * est.guess
		guess[ paste(lav1$lhs) ] <- lav1$ustart	

		names(guess) <- names(est.guess) <- NULL
		
		res$guess <- guess
		res$est.guess <- est.guess
		res$method <- "tam.mml.3pl"
	
						}
							
	return(res)
		}
