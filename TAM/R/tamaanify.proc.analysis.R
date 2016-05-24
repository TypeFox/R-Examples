



#########################################################
# process analysis
tamaanify.proc.analysis <- function( res ){ 	
	tam1 <- res$tammodel.dfr
	ind1 <- which( paste(tam1$syn) == "ANALYSIS:" )
	res$ANALYSIS <- ""	
	if ( length(ind1) > 0 ){
		index1 <- tam1$part_begin[ ind1 ]
		dfr <- paste( tam1[ which( tam1$part_begin == index1 )[-1] , "syn" ] )
		res$ANALYSIS <- dfr	
		ANALYSIS.list <- list()	
		if( length( grep("TYPE=LCA",dfr) ) > 0 ){
			ANALYSIS.list$type <- "LCA"
					}
		if( length( grep("TYPE=MIXTURE",dfr) ) > 0 ){
			ANALYSIS.list$type <- "MIXTURE"
					}		
		# located latent class analysis			
		if( length( grep("TYPE=LOCLCA",dfr) ) > 0 ){
			ANALYSIS.list$type <- "LOCLCA"
					}		
		# trait		
		if( length( grep("TYPE=TRAIT",dfr) ) > 0 ){
			ANALYSIS.list$type <- "TRAIT"
					}						
		# ordered latent class analysis			
		if( length( grep("TYPE=OLCA",dfr) ) > 0 ){
			ANALYSIS.list$type <- "OLCA"
					}
		ind <- grep("NCLASSES" , dfr )			
		if ( length(ind) > 0 ){
			m1 <- gsub( "NCLASSES(" , "" , dfr[ind] , fixed=TRUE )
			m1 <- as.numeric( gsub(")" , "" , m1 , fixed=TRUE ) )
			ANALYSIS.list$NCLASSES <- m1
							}
		ind <- grep("NSTARTS" , dfr , fixed=TRUE)			
		if ( length(ind) > 0 ){
			m1 <- gsub( "NSTARTS(" , "" , dfr[ind] , fixed=TRUE )
			m1 <- gsub(")" , "" , m1 , fixed=TRUE )
			m1 <- as.numeric( unlist( strsplit( m1 , split= "," , fixed=TRUE) ) )
			ANALYSIS.list$NSTARTS <- m1
							}				
					
		res$ANALYSIS.list <- ANALYSIS.list
		} else {
		
		res$ANALYSIS.list$type <- "TRAIT"
			}
		res$skillspace <- "normal"	
		return(res)		
		}
##########################################################
