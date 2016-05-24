

					
					

.generate.interactions <- function(X , facets , formulaA , mm ){	
	d1 <- d0 <- X	
	h1 <- sapply( colnames(d1) , FUN = function(vv){
				length(grep( vv , paste(formulaA) )) } )
	h1 <- colnames(d1)[ h1 == 0 ]
	d0 <- d0[ , ! ( colnames(d1) %in% h1 ) , drop=FALSE]
	M2 <- stats::model.matrix( #object= 
			formulaA , data= d1 , 
			contrasts.arg = lapply( d0 , contrasts, contrasts=FALSE) )
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
					varsc <- na.omit( varsc)
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
######################
# rename item names
.rename.items <- function( matr , itemren , cols=TRUE ){
	rM <- rownames(matr)
	cM <- colnames(matr)
	I <- nrow(itemren)

	vers <- FALSE
	vers <- TRUE

	#-----------------------
	# rows
	if (vers){
		rM0 <- rM
		itemren2 <- paste0(itemren[,2] , "-")
		nc2 <- nchar(itemren2)
		N1 <- min(nc2)
		N2 <- max(nc2)
		for (nn in N1:N2){
			i1 <- match( substring( rM0 , 1 , nn ) , itemren2 )
			h1 <- paste0( itemren[ i1 ,1] , "-"  , substring( rM0 , nn+1 , nchar(rM0) ) )
			i2 <- ! is.na(i1)
			rM0[ i2 ] <- h1[ i2]
						  }
		rM <- rM0
		if ( cols){
			cM0 <- cM	
			ind <- match(  cM0 , itemren[,2])		
			ind <- stats::na.omit(ind)
			cM0[ ind ] <- paste(itemren[,1])	
			itemren2 <- paste0(itemren[,2] , ":")
			nc2 <- nchar(itemren2)
			N1 <- min(nc2)
			N2 <- max(nc2)
			for (nn in N1:N2){
				i1 <- match( substring( cM0 , 1 , nn ) , itemren2 )
				h1 <- paste0( itemren[ i1 ,1] , ":"  , substring( cM0 , nn+1 , nchar(cM0) ) )
				i2 <- ! is.na(i1)
				cM0[ i2 ] <- h1[ i2]
							  }
			cM <- cM0
					}
					
				}

if (!vers){	
	for ( ii in 1:I){
		rM <- gsub( paste0( itemren[ii,2] , "-") , paste0( itemren[ii,1] , "-") , rM )
		if (cols){
			cM <- gsub( paste0( itemren[ii,2] , ":") , paste0( itemren[ii,1] , ":") , cM )	
			cM[ cM == itemren[ii,2] ] <- paste(itemren[ii,1])
				}
					}					
			}
			
	rM -> rownames(matr)

	if ( cols){ cM -> colnames(matr) }
	return(matr)
		}
#############################################################		
.rename.items2 <- function( vec , itemren ){
	cM <- vec
	I <- nrow(itemren)

vers <- TRUE	
#vers <- FALSE
if (vers){
v0 <- Sys.time()
	cM0 <- cM
	# relabel items
	ind <- match(  cM , itemren[,2])
	ind <- na.omit(ind)
	cM0[ ind ] <- paste(itemren[,1])
	itemren2 <- paste0(itemren[,2] , ":")
	nc2 <- nchar(itemren2)
	N1 <- min(nc2)
	N2 <- max(nc2)
	for (nn in N1:N2){
		i1 <- match( substring( cM0 , 1 , nn ) , itemren2 )
		h1 <- paste0( itemren[ i1 ,1] , ":"  , substring( cM0 , nn+1 , nchar(cM0) ) )
		i2 <- ! is.na(i1)
		cM0[ i2 ] <- h1[ i2]
					  }
	cM <- cM0	
	}
if (!vers){	
	for ( ii in 1:I){
			cM <- gsub( paste0( itemren[ii,2] , ":") , paste0( itemren[ii,1] , ":") , cM )	
			cM[ cM == itemren[ii,2] ] <- paste(itemren[ii,1])
				}
			}
	return(cM)
		}
#############################################################
.rename.items3 <- function( matr , facet.list , I , cols=TRUE  ){
### check for equalities in rM and cM in all entries!!!!
	rM <- rownames(matr)
	rMsplit <- strsplit( paste(rM) , split="-" )	
	RR <- length(rMsplit)
	FF <- length(facet.list)
	for (rr in 1:RR){
		rr1 <- rMsplit[[rr]]
		if (FF > 0 ){
		for (ff in 1:FF){ # begin ff
#		for (ff in seq(1,FF,1) ){
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)
				for (ii in 1:I ){ 
					rr1[ rr1 == itemren[ii,2] ] <- paste(itemren[ii,1])
					rMsplit[[rr]] <- rr1
								}
							} # end ff
							} # end if FF
					}
	rM <- unlist( lapply( rMsplit , FUN = function(ll){ paste( ll , collapse="-") } )	)
	rownames(matr) <- rM	
	#****************************************
	if ( cols){
		cM <- colnames(matr)
		cMsplit <- strsplit( paste(cM) , split=":" )	
		RR <- length(cMsplit)
		FF <- length(facet.list)
		for (rr in 1:RR){
			rr1 <- cMsplit[[rr]]
			if (FF>0){
			for (ff in 1:FF){ # begin ff
				itemren <- facet.list[[ff]]
				I <- nrow(itemren)
					for (ii in 1:I ){ 
						rr1[ rr1 == itemren[ii,2] ] <- paste(itemren[ii,1])
						cMsplit[[rr]] <- rr1
									}
								} # end ff
							}
						}
		cM <- unlist( lapply( cMsplit , FUN = function(ll){ paste( ll , collapse=":") } )	)
		colnames(matr) <- cM	
			}
	return(matr)
		}
#############################################################
.rename.items2a <- function( vec , facet.list , I ){
### check for equalities!!!
	cM <- vec
	FF <- length(facet.list)
	rM <- cM
v0 <- Sys.time()	
    if ( ! is.null(rM) ){ 
		rMsplit <- strsplit( paste(rM) , split="-" )	
		RR <- length(rMsplit)
		FF <- length(facet.list)
		NRM <- max( unlist( lapply( rMsplit , FUN=function(ll){ length(ll) } ) ))
		rMsplit0 <- rMsplit
		rMsplit0 <-  matrix( unlist( rMsplit0 )  , ncol=NRM , byrow=TRUE )
# cat(" *** ren2a split " ) ; v1 <- Sys.time() ; print(v1-v0) ; v0 <- v1 
		if (FF>0){
		  for (ff in 1:FF){
			itemren <- facet.list[[ff]]
			for (nn in 1:NRM){
			# nn <- 3
			rm_nn <- rMsplit0[,nn]
			ind1 <- match( rm_nn , itemren[,2] )
			ind1 <- stats::na.omit(ind1)
			h1 <- paste(itemren[ ind1 , 1] )
			if ( length(h1) > 0 ){
				rMsplit0[,nn] <- h1 
						}
				}
					}
# cat(" *** loop facets " ) ; v1 <- Sys.time() ; print(v1-v0) ; v0 <- v1 
			cM0 <- apply( rMsplit0 , 1 , FUN = function(ll){ paste0( ll , collapse="-") } )
			cM <- cM0
# cat(" *** apply " ) ; v1 <- Sys.time() ; print(v1-v0) ; v0 <- v1 			
				}
		}
			
			
	return(cM)
		}
#############################################################
.rename.items2aa <- function( vec , facet.list , I ){
### check for equalities!!!
	cM <- vec
	FF <- length(facet.list)
	rM <- cM

    if ( ! is.null(rM) ){ 
		rMsplit <- strsplit( paste(rM) , split="-" )	
		RR <- length(rMsplit)
		FF <- length(facet.list)
		NRM <- max( unlist( lapply( rMsplit , FUN=function(ll){ length(ll) } ) ))
		rMsplit0 <- rMsplit
		rMsplit0 <-  matrix( unlist( rMsplit0 )  , ncol=NRM , byrow=TRUE )

		if (FF>0){
		  for (ff in 1:FF){
			itemren <- facet.list[[ff]]
			for (nn in 1:NRM){
			# nn <- 3
			rm_nn <- rMsplit0[,nn]
			ind1 <- match( rm_nn , itemren[,2] )
			ind1 <- stats::na.omit(ind1)
			h1 <- paste(itemren[ ind1 , 1] )
			if ( length(h1) > 0 ){
				rMsplit0[,nn] <- h1 
						}
				}
					}

			cM0 <- apply( rMsplit0 , 1 , FUN = function(ll){ paste0( ll , collapse="-") } )
			cM <- cM0
				}
		}
	return(cM)
		}
##################################################################

			
		
#############################################################		
#############################################################
.rename.items3a <- function( matr , facet.list , I , cols=TRUE ,
			xsi.table ){		
### check for equalities in rM and cM in all entries!!!!
	rM <- rownames(matr)
	rMsplit <- strsplit( paste(rM) , split="-" )
	RR <- length(rMsplit)
	rMM <- matrix( unlist(rMsplit) , nrow=RR , byrow=TRUE)
	rMM.ncol <- ncol(rMM)
	FF <- length(facet.list)	
	
	
	rMM1 <- rMM	
	NRM <- ncol(rMM1)
	if (FF>0){
	for (ff in 1:FF){
		# ff <- 1  # facet ff
		itemren <- facet.list[[ff]]
		for (nn in 1:NRM){
		   ind1 <- match( rMM1[,nn] , itemren[,2] )
		   ind1 <- stats::na.omit(ind1)
		   h1 <- paste(itemren[ ind1 , 1] )
		   if ( length(h1) > 0 ){
				rMM1[,nn] <- h1 
						}
				}
				}		
			}		
	rMM1 <- apply( rMM1 , 1 , FUN = function(ll){ paste0( ll , collapse="-") } )
	rM <- rMM1		
		
		
	
	if ( cols){
		rM <- colnames(matr)		
		rMsplit <- unlist( strsplit( paste(rM) , split=":" ) )
	    xsi.table <- xsi.table[xsi.table$constraint==0,]
		XT <- nrow(xsi.table)
		F0 <- max(xsi.table$facet.order)				
		index <- sapply( 1:XT , FUN = function(xx){
				m1 <- cbind( xx , 1:xsi.table[xx,"facet.order"] )
				matrix( t(m1) , ncol=1 , byrow=FALSE)
								} )					
		index <- matrix( unlist(index) , ncol=2 , byrow=T)
		rMMsub <- matrix("" , nrow=XT, ncol=F0) 
		rMMsub[ index ] <- rMsplit
		FF <- length(facet.list)	
		if (FF>0){ 
		for (ff in 1:FF){ # ff <- 1
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)		
			for (ii in 1:I){ # ii <- 1
			for (kk in 1:F0){# kk <- 3
				rMMsub[ rMMsub[,kk] == itemren[ii,2] , kk ] <- paste(itemren[ii,1])
								}
							}
						}
					}	# end if FF>0
		rM <- unlist( sapply( 1:XT , FUN = function(kk){
			paste( rMMsub[ kk , seq(1 , xsi.table$facet.order[kk] ) ] , collapse=":" ) } )
						)
		colnames(matr) <- rM					
			}
	return(matr)
		}
################################################################################		
################################################################################
.rename.items2b <- function( vec , facet.list , I , xsi.table , sel1=0 ){
### check for equalities!!!
#	cM <- vec	
	rM <- vec
	if ( ! is.null(rM)){
		rMsplit <- unlist( strsplit( rM , split=":" ) )
		if (sel1==1){ xsi.table <- xsi.table[xsi.table$constraint==0,] }
		if (sel1==2){ xsi.table <- xsi.table[xsi.table$constraint==1,] }		
		XT <- nrow(xsi.table)
		F0 <- max(xsi.table$facet.order)
		index <- sapply( 1:XT , FUN = function(xx){
				m1 <- cbind( xx , 1:xsi.table[xx,"facet.order"] )
				matrix( t(m1) , ncol=1 , byrow=FALSE)
								} )					
		index <- matrix( unlist(index) , ncol=2 , byrow=T)
		rMMsub <- matrix("" , nrow=XT, ncol=F0) 
		rMMsub[ index ] <- rMsplit
		FF <- length(facet.list)	
		if (FF>0){ 
		for (ff in 1:FF){ # ff <- 1
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)		
			for (ii in 1:I){ # ii <- 1
			for (kk in 1:F0){# kk <- 3
				rMMsub[ rMMsub[,kk] == itemren[ii,2] , kk ] <- paste(itemren[ii,1])
								}
							}
						}
					}	
		rM <- unlist( sapply( 1:XT , FUN = function(kk){
			paste( rMMsub[ kk , seq(1 , xsi.table$facet.order[kk] ) ] , collapse=":" ) } )
						)
				}
	
	return(rM)
		}
		
		
