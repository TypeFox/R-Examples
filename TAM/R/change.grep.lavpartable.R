

		
################################################################
# change lavaan parameter label
change.grep.lavpartable <- function( lavpartable ){
    labels1 <- lavpartable$label
    labels1 <- sort( unique( labels1 ) )	
    labels1 <- grep( "__" , labels1 , value=TRUE  )
	labels1 <- labels1[ substring( labels1 , nchar(labels1)-2 , nchar(labels1) ) != "__" ]	
    LL <- length(labels1)
    if (LL > 0 ){
        for (ll in 1:LL){
            # ll <- 1
            label.group <- labels1[ll]	
            lav.changed <- extend.label.group( label.group )
            ind.ll <- which( lavpartable$label %in% label.group )
			if ( length(ind.ll) == length(lav.changed) ){ 
				lavpartable[ ind.ll , "label" ] <- lav.changed 
									}
                        }
                }
	res <- list("lavpartable" = lavpartable , "changed" = LL > 0 )				
    return(res)
            }			
