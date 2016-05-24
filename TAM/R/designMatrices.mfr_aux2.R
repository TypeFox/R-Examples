

#######################################################
# generate interactions										
.generate.interactions2 <- function(X , facets , formulaA , mm ){	
	d1 <- d0 <- X	
	h1 <- sapply( colnames(d1) , FUN = function(vv){
				length(grep( vv , paste(formulaA) )) } )
	h1 <- colnames(d1)[ h1 == 0 ]
	d0 <- d0[ , ! ( colnames(d1) %in% h1 ) , drop=FALSE]
	M2 <- stats::model.matrix( #object= 
         formulaA , data= d1 , 
			contrasts.arg = lapply( d0 , stats::contrasts, contrasts=FALSE) )
	h2 <- colnames(M2)
	h1 <- colnames(mm)
	# extract facets
	xsi.table <- data.frame( "parameter" = h2 )
	xsi.split <- sapply( xsi.table$parameter , FUN = function(ll){ 
		l1 <- as.vector( unlist( strsplit( paste(ll) , split=":" ) ) )
		v1 <- l1
		for (ii in 1:length(l1) ){
			for (cc in colnames(X) ){ 
				kk <- grep( cc , l1[ii] )
				if (length(kk)>0){ v1[ii] <- cc }
								}
						} 
		v1 <- paste0( v1 , collapse=":" )
		return(v1)		
					} )
	xsi.table$facet <- unlist(xsi.split)
	xsi.table$facet.order <- sapply( xsi.table$parameter , FUN = function(ll){ 
		length( as.vector( unlist( strsplit( paste(ll) , split=":" ) ) ) ) } )
	xsi.table$constraint <- 1 - 1*(xsi.table$parameter %in% h1)
	xsi.table$facet.index <- match( xsi.table$facet , unique( xsi.table$facet ) )
#	xsi.table$orig.index <- seq(1,nrow(xsi.table))
#	xsi.table[ order( paste( xsi.table$facet.index+100 , xsi.table$parameter ) ) , ]

	facets.unique <- unique( xsi.table$facet )
	b1 <- xsi.table[ xsi.table$constraint == 1 , "parameter" ]
	c1 <- xsi.table[ xsi.table$constraint == 0 , "parameter" ]	
    xsi.constraints <- matrix( NA  , nrow=length(b1) , ncol=length(c1) )
	rownames(xsi.constraints) <- paste(b1)
	colnames(xsi.constraints) <- paste(c1)
# Revalpr("xsi.constraints")
# Revalpr("b1")
# stop()	
	
	
# b1 <- b1[3]
	############################
	# loop over terms
    for (bb in b1 ){
		#bb <- b1[9]
		v1 <- 0
		mult <- 1
		xsi.table.bb <- xsi.table[ xsi.table$parameter == bb , ]
		x0 <- x1 <- xsi.table[ xsi.table$facet %in% xsi.table.bb$facet , ]
		if ( xsi.table.bb$facet.order==1){ 
			xsi.constraints[paste(bb),] <- 0		
			xsi.constraints[ paste(bb) , paste( x1[ x1$constraint == 0 , "parameter" ] ) ] <- -1
										}
		if ( xsi.table.bb$facet=="item:step"){ 
			v1 <- 1
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					ll[1] == s2[1] } ) )
			x1 <- x1[ g1 , ]	
			mult <- 1
# cat("......",bb,"......\n")
# print(x1)			
			varsc <- paste( x1[ x1$constraint == 0 , "parameter" ] )
			if ( length(varsc) == 0){
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						ll[2] == s2[2] } ) )
				x1 <- x0[ g1 , ]	
				varsc <- paste(x1[ x1$constraint == 0 , "parameter" ])
				mult <- 1
				if ( length(varsc) == 0){
					varsc <- x1[ , "parameter" ]
					varsc <- setdiff( varsc , paste(bb) )
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					mult <- -1
									}				
						}
			
			xsi.constraints[ paste(bb) , varsc ] <- -1*mult
										}
		##################
		### order 2
		if ( xsi.table.bb$facet.order==2 & v1 ==0 ){ 							
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )	
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					ll[2] == s2[2] } ) )
			x1 <- x1[ g1 , ]	
			varsc <- x1[ x1$constraint == 0 , "parameter" ]
			if ( length(varsc) == 0){
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						ll[1] == s2[1] } ) )
				x1 <- x0[ g1 , ]	
				varsc <- x1[ x1$constraint == 0 , "parameter" ]	
				if ( length(varsc) == 0){
					varsc <- x1[ , "parameter" ]
					varsc <- setdiff( varsc , paste(bb) )
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					mult <- -1
									}
								}							
			xsi.constraints[ paste(bb) , paste( varsc ) ] <- -1	* mult									
										}
		#########################
		### order 3
		if ( xsi.table.bb$facet.order==3 & v1 ==0 ){ 			
			mult <- 1	
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )	
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					( ll[2] == s2[2] ) & (ll[3] == s2[3] )  } ) )
			x1 <- x1[ g1 , ]
			varsc <- x1[ x1$constraint == 0 , "parameter" ]
			if ( length(varsc) == 0 ){			
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						( ll[1] == s2[1] ) & (ll[3]==s2[3]) } ) )
				x1 <- x0[ g1 , ]	
				varsc <- x1[ x1$constraint == 0 , "parameter" ]		
				
				if ( length(varsc) == 0 ){			
					varsc <- x1[ , "parameter" ]			
					varsc <- setdiff( varsc , paste(bb) )	
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					varsc <- stats::na.omit( varsc)
					mult <- -1			
						if ( length(varsc) == 0 ){
							g1 <- unlist( lapply( s20 , FUN = function(ll){
									( ll[1] == s2[1] ) & (ll[2]==s2[2]) } ) )
							x1 <- x0[ g1 , ]	
							varsc <- x1[ x1$constraint == 0 , "parameter" ]		
							mult <- 1	
							
							if ( length(varsc) == 0 ){ 	
#							varsc <- setdiff( varsc , paste(bb) )		
								varsc <- x1[ , "parameter" ]
   							    varsc <- setdiff( varsc , paste(bb) )	
						        h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)					
								varsc <- names(h1)[  h1 != 0 ]
								varsc <- stats::na.omit( varsc)
								mult <- -1					
													}
						}
										}																
							}		
						
			if ( length(varsc) > 0 ){
				xsi.constraints[ paste(bb) , paste( varsc ) ] <- -1	* mult	
					} else {
				xsi.constraints[ paste(bb) , ] <- NA	
						}
	
		
							}
						}
		xsi.constraints[ rowSums( abs(xsi.constraints) ) == 0 , ] <- NA
		res <- list( "xsi.constraints" = xsi.constraints , "xsi.table" = xsi.table )
#print(res) ;  stop("here")		
		
		return(res)
}						
