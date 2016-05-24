
###################################################
# define R function (method)
tamaanify.define.method <- function(res , tam.method ){

	itemtypes <- paste0( res$items$itemtype )
	l1 <- strsplit( itemtypes , split="," , fixed=TRUE )
	itemtypes <- unlist( lapply( l1 , FUN = function(ll){ ll[1] } ) )

	m1 <- mean(  itemtypes %in% c("Rasch" , "PCM" ) )
	
	if ( m1 == 1 ){  
			res$method <- "tam.mml" 
					}
	if ( ! is.null(tam.method) ){
		res$method <- tam.method 
					}				
	al <- res$ANALYSIS.list$type
	if ( al %in% c("LCA", "LOCLCA","MIXTURE" , "OLCA") ){
		res$method <- "tam.mml.3pl"
						}	
	#--------------------------------
    #**** choose generalized partial credit model
	if ( res$method %in% c("tam.mml.2pl" ) ){	
         items <- res$items
		 m1 <- mean( items$itemtype %in% c( "GPCM" ) )
		 if (m1==1){
		     res$irtmodel <- "GPCM"
					}
					}
					
	return(res)
			}
###################################################			

