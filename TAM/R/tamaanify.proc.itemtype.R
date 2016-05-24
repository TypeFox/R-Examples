


############################################################
# process item type
tamaanify.proc.itemtype <- function( res ){
    tam1 <- res$tammodel.dfr	
	ind1 <- which( paste(tam1$syn) == "ITEMTYPE:" )
	
	if ( length(ind1) > 0 ){
		index1 <- tam1$part_begin[ ind1 ]
		ind2 <- which( tam1$part_begin == index1 )[-1]
		if ( length(ind2) > 0 ){
			m2 <- paste( tam1[ ind2 , "syn" ] )
			res$ITEMTYPE <- m2
			#*** grep item types	
			m2 <- data.frame( "index" = seq(1,length(m2) ) , "syn" = m2 )
			m3 <- strsplit( paste( m2$syn) , split="(" , fixed=TRUE )
			m2$items <- unlist( lapply( m3 , FUN = function(ll){ ll[1] } ) )
			m2$itemtype <- unlist( lapply( m3 , FUN = function(ll){ ll[2] } ) )
			m2$itemtype <- gsub( ")" , "" , paste(m2$itemtype) , fixed=TRUE )
			N2 <- nrow(m2)
			items0 <- res$items
			items00 <- paste( items0$item )
			for (ii in 1:N2){
				# ii <- 1
				h1 <- unlist(strsplit( paste(m2$items[ii]) , split="__" ))
				if ( h1[1] == "ALL" ){
				
					h1 <- c( items00[1] , items00[ length(items00) ] )
									}
									
				if ( length(h1) == 1 ){
					h1 <- c( h1 , h1 )
										}
				ind1 <- which( items00 == h1[1] )
				ind2 <- which( items00 == h1[2] )
				items0[ seq(ind1,ind2) , "itemtype"] <- paste(m2$itemtype[ii])
									}
				res$items <- items0						
							}											
					}

     return(res)
		}
############################################################		
