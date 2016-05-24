


##########################################################
# grep list of linear equations
tamaanify.grep.linequations <- function( syn2 ){		
		# split syntax in matrix format				
		syn3 <- strsplit( paste0(syn2) , split="==" , fixed=TRUE )
		syn3a <- unlist(lapply( syn3 , FUN = function(ll){ ll[1] } ))
		syn3b <- unlist(lapply( syn3 , FUN = function(ll){ ll[2] } ))
		L1 <- length(syn3a)
		dfr <- NULL
		for (ii in 1:L1){
			# ii <- 2
			ll1 <- unlist( strsplit( syn3b[ii] , split="+" , fixed=TRUE ) )
			dfr2 <- data.frame("lhsparm" = syn3a[ii] , "rhsparm" = syn3b[ii] ,
					"terms" = ll1 )
			dfr <- rbind( dfr , dfr2 )
						}

		m2 <- strsplit( paste(dfr$terms)	, split="*" , fixed=TRUE )
		# numerical
		dfr$fac <- as.numeric( unlist( lapply( m2 , FUN = function(uu){
					m1 <- ifelse( length(uu) == 1 , 1 , uu[1] )
					m1 <- gsub( "(" , "" , m1 , fixed=TRUE )
					m1 <- gsub( ")" , "" , m1 , fixed=TRUE )
					m1
							} ) ) )
		# parameter
		dfr$parm <- unlist( lapply( m2 , FUN = function(uu){
					ifelse( length(uu) == 1 , uu[1] , uu[2] )
							} ) )
		rownames(dfr) <- NULL							
		return(dfr)
		}
###############################################################		

