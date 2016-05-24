


############################################################
# process item characteristics
tamaanify.proc.items <- function( res , resp){
    lavpartable <- res$lavpartable

	ind <- which( lavpartable$op == "=~" )
	items0 <- unique( lavpartable[ind,"rhs"] )
	
	maxK <- apply( resp[,items0] , 2 , max , na.rm=TRUE )
	items <- data.frame( "item" = names(maxK) , "ncat" = maxK+1 )
	items <- items[ items$item %in% items0 , ]
	items$itemtype <- ifelse( items$ncat == 2 , "2PL" , "GPCM" )
	if( res$ANALYSIS.list$type == "LOCLCA" ){
		items$itemtype <- ifelse( items$ncat == 2 , "Rasch" , "PCM" )
								}
	items$itemtype <- gsub("RASCH" , "Rasch" , items$itemtype )
	res$items <- items
	# res$maxcat <- max(maxK)
	res$maxcat <- max( maxK[ names(maxK) %in% items0 ] )
	
	
	#********************************
	# items with guessing parameters
	
	l2 <- lavpartable[ lavpartable$op == "?=" , ]
	if ( nrow(l2) > 0 ){
		items <- res$items
		items[ items$item %in% l2$lhs , "itemtype" ] <- "3PL"
		res$items <- items
						}
								
	# add all thresholds in lavaan parameter table
	lavpartable$fullsyn <- paste0( lavpartable$lhs , lavpartable$op , lavpartable$rhs )	
	lavpartable2 <- lavpartable		
	lav0 <- lavpartable[1,]

	#****** search for items with missing thresholds in parameter table
	for ( hh in 1:(max(maxK)) ){		
		# hh <- 1			
		types <- paste0( items$item , "|t" , hh )
		ind <- which( ! ( types %in% lavpartable$fullsyn ) )
		LI <- length(ind)
		if (LI > 0 ){
			lav1 <- lav0[ rep(1,LI) , ]
			lav1$fullsyn <- types[ind]
			lav1$label <- ""
			lav1$op <- "|"
			lav1$rhs <- paste0( "t" , hh )
			lav1$lhs <- items$item[ind]
			lav1$id <- 999
			lav1$free <- 999			
			item.sel <- items$item[ind]
			m1 <- paste(items$item[ ( items$ncat - 1 ) < hh ] )
			if ( length(m1) > 0 ){
				lav1[ lav1$lhs %in% m1 , "free" ] <- 0
				lav1[ lav1$lhs %in% m1 , "ustart" ] <- 999
				lav1[ lav1$lhs %in% m1 , "user" ] <- -99
								}
			lavpartable2 <- rbind( lavpartable2 , lav1 )
						}	
			}
	
	rownames(lavpartable2) <- NULL
	
	#********************
	# define new parameter labels in parameter table

	#*** loadings
	ind <- which( ( paste(lavpartable2$op) == "=~" ) & ( paste( lavpartable2$label ) == "" ) )
	if ( length(ind) > 0 ){
		lavpartable2[ ind , "label" ] <- paste0( lavpartable2$rhs[ind] , "_" , 
					lavpartable2$lhs[ind]  , "_load" )
							}
	#*** thresholds
	for (hh in 1:( max(maxK) ) ){	
		ind <- which( ( paste(lavpartable2$rhs) == paste0("t" , hh) ) & ( paste( lavpartable2$label ) == "" ) &
						( paste(lavpartable2$op) == "|" )  )					
		if ( length(ind) > 0 ){
			lavpartable2[ ind , "label" ] <- paste0( lavpartable2$lhs[ind] , "_" , "Cat" , hh )
								}
							}				
	res$lavpartable <- lavpartable2
	
	return(res)	
		}
