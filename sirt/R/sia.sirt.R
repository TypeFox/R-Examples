
###################################################
# statistical implicative analysis
# simplified algorithm
sia.sirt <- function(dat , significance=.85 ){

	dat1 <- dat
	# order items
	dat1 <- dat1[ , order( colMeans(dat1 , na.rm=TRUE ) ) ]	
	p <- colMeans(dat1 , na.rm=TRUE )
	
	
	# handle missing responses
	dat1.resp <- 1*(1-is.na(dat1) )
	dat1[ is.na(dat1) ] <- 0
	I <- ncol(dat1)
	
	# total sample size
	ntot <- crossprod( dat1.resp )
	# n10
	n10 <- crossprod( dat1.resp * ( dat1 == 1 ) , dat1.resp * ( dat1 == 0 ) )
	# n11
	n11 <- crossprod( dat1.resp * ( dat1 == 1 ) , dat1.resp * ( dat1 == 1 ) )
	# calculate some probabilities
	p10 <- n10 / ntot       # X_i=1, X_j = 0
	p11 <- n11 / ntot
	# item p values in every cell
	p1 <- diag( n11 ) / diag( ntot )
	p1M <- outer( p1 , p1 )
	p0 <- 1-p1
	p0M <- 1 - p1M
	# probability under independence
	p10_ind <- outer( p1 , p0 )
	# define effect size
	impl_int_es <- ( p10 - p10_ind  ) / sqrt( p10_ind )
	diag(impl_int_es) <- NA
	# implicative intensity t value
	impl_int_t <- sqrt( ntot ) * impl_int_es
	# significance value: implicative intensity
	impl_int <- 1 - stats::pnorm( impl_int_t )
	# look at significant implications
	impl_significance <- 1 * ( impl_int > significance )
	conf_loev <- p11 / outer( p1 , rep(1,I) )
	diag( conf_loev) <- colMeans( dat1 , na.rm=TRUE )

	#***
	# compute arrows: remove symmetric and transitive paths
	I1 <- impl_significance 
	# remove symmetric implications
#	I1[ ( I1 == t(I1) ) & ( I1 == 1 ) ] <- 0
	for (ii in 1:(I-1) ){
		for (jj in (ii+1):I){
			   l2 <- I1[ii,jj] + I1[jj,ii]
			   if (l2 == 2 ){
					if (impl_int[ii,jj] > impl_int[jj,ii] ){
							I1[jj,ii] <- 0 } else { I1[ii,jj] <- 0 }
			   
							} # end l2 = 2
			   
							}
						}


# print(I1)
# stop()						
						
	# remove transitive relations
	I1 <- .sia.remove.transitive(I1)	
	
	# calculate potencies to look for connected graphs
	I1.pot <- I1
	
	# calculate matrix potence
	IS <- IS0 <- I1
	pp <- 2
	dev <- 1
	while ( dev > .5 ){
		IS.old <- IS
		G1 <-  1 * ( IS0 %*% I1 > 0 )
		IS0 <- G1
		I1.pot[ G1 > 0 ] <- pp
		IS <- 1 * ( IS + G1  > 0 )
		pp <- pp + 1
		dev <- sum ( abs( IS - IS.old ) )
#		cat( "Iteration " , pp , "dev " , dev , "\n")
					}
	
	# descriptives
	desc <- list( "nodes" = I  )
	# unconnected nodes	
	ind.uc <- intersect( which( rowSums(I1) == 0 ) , which( colSums(I1) == 0 ) )
	UC <- length( ind.uc )
	vars.uc <- rownames(I1)[ ind.uc ]
	if (  UC > 0 ){
		desc$nodes.unconnected <- UC
					}
	desc$edges <- sum(I1)					
	desc$sign.implications <- sum( impl_significance , na.rm=TRUE) 	
	# descriptives on item level
	desc.item <- data.frame( "item" = colnames(dat1) , "p" = p )						
	desc.item$level <- g1 <- apply( I1.pot , 1 , max )
# cat("\nc100")				
	#****
	# recalculate levels
	g0 <- g1	
    g0 <- sort( g0 , decreasing = TRUE )
    g1 <- g0	
    for (ii in 1:(I-1) ){
		# ii <- 1
		vii <- names(g1)[ii]
		levii <- g1[ vii ]
		hzu <- which( I1[vii,] > 0 )
		if ( length( hzu) > 0 ){
			g1[ colnames(I1)[hzu] ] <- levii - 1
			g1 <- sort( g1 , decreasing = TRUE )
							}
					}
	desc.item[ names(g1) , "level" ]  <- g1
	if (UC > 0 ){ 
#		desc.item$level[ind.uc] <- NA 
#		g1[ind.uc] <- -1
		desc.item[vars.uc,"level"] <- NA 
		g1[vars.uc] <- -1
			}
	# from objects and to objects
	desc.item$from <- colSums(I1)
	desc.item$to <- rowSums(I1)
	desc.item <- desc.item[ order( desc.item$level , decreasing = TRUE ) , ]					
#	rownames(desc.item) <- NULL	
	
	# create graph matrix
	dfr <- NULL
	dfr2 <- NULL
    g1 <- g1[ colnames(I1) ]
	CI1 <- paste0( colnames(I1) , "\n(L" , g1 , ")")
	CI2 <- colnames(I1)
	for (ii in 1:I){
		# ii <- 1	
		v1 <- which( I1[ii, ] > 0 )
		if ( length(v1) > 0 ){
			dfr.ii <- cbind( CI1[ii] , CI1[ v1] )
			dfr <- rbind( dfr , dfr.ii )
			dfr2.ii <- cbind( CI2[ii] , CI2[ v1] )
			dfr2 <- rbind( dfr2 , dfr2.ii )			
						}
				}
	if (UC > 0 ){
			v1 <- ind.uc
			dfr.ii <- cbind( CI1[v1] , CI1[ v1] )
			dfr <- rbind( dfr , dfr.ii )
			dfr2.ii <- cbind( CI2[v1] , CI2[ v1] )
			dfr2 <- rbind( dfr2 , dfr2.ii )			
					}
					
	# create igraph object
	igraph.obj <-  igraph::graph.edgelist(dfr)
	# create Rgraphviz object
	
	
	
	#*****
	# OUTPUT
	res <- list( "adj.matrix" = I1 , 
	    "adj.pot" = I1.pot , 
		"adj.matrix.trans" = IS , 
		"desc" = desc , "desc.item" = desc.item ,
		"impl.int" = impl_int , "impl.t" = impl_int_t , 
		"impl.significance" = impl_significance , 
		"conf.loev" = conf_loev , 
		"graph.matr" = dfr2 , 
		"graph.edges" = unique( c( dfr2[,1] , dfr2[,2]  ) ) ,
		"igraph.matr" = dfr , "igraph.obj" = igraph.obj)
	return(res)
		
		}
##################################################################


####################################################################
# remove transitive relations
.sia.remove.transitive <- function(I1){
	I <- ncol(I1)	
	diag(I1) <- 0	
	BB <- 1
	IS <- I1
	iter <- 0
	while (BB > .0001){
		I0 <- IS
		iter <- iter + 1
		for (ii in 1:I ){
			for (jj in 1:I){
				if (ii!= jj ){
					if ( sum( IS[ii,] * IS[,jj] ) > 0 ){ 
							I1[ii,jj] <- 0
									}
						}
					}
				}
		  IS <- IS + IS %*% I1
		  IS <- 1*(IS > 0 )
		  BB <- sum( abs( I0 - IS ) )
#		  BB <- 0
			}
			
	return(I1)
		}
####################################################################