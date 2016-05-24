
# Examine relationship between found phrases
#


#' Make phrase appearance matrix from textreg result.
#'
#' Make matrix of which phrases appear in which of the positively
#' marked documents.
#'
#' Very similar to phrase.matrix, except this looks only at positively 
#' marked documents and just returns 1 or 0 on whether any document has
#' a phrase, rather than giving counts.
#' This is used by the clustering vizualizations and make.similarity.matrix.
#'
#' @seealso make.similarity.matrix
#' @seealso phrase.matrix
#'
#' @param result  An textreg.result object.
#' @return A $n X p$ matrix for $n$ documents and $p$ phrases in the result object.  Each entry is a 0/1 value 
#'    indicating presence of the given phrase in the given document.
#' @family Phrase Vizualization
#' @export
make.appearance.matrix = function( result ) {
	mat = phrase.matrix(result)
	mat = mat[ result$labeling==1, ]
	#head(mat)
	
	stopifnot( all( result$model$ngram == colnames(mat) ) )
	# normalize by Z
	#mat = sweep( mat, 2, result$model$Z, FUN="/" )
	
	# turn matrix into appearance of phrase in documents.
	mat[ mat > 0 ] = 1
	
	# drop intercept
	#mat = mat[,-1]

	# if single column, make sure we are real matrix
	if ( !is.matrix( mat ) ) {
		mat = matrix( mat, ncol=1, nrow=length(mat) )
	}
	
	mat
}
	
	
#' Calculate similarity matrix for set of phrases.
#'
#' First get phrase appearance pattern on positive labeling 
#' (if not directly passed)
#' and then calculate similarity matrix of how they are similar
#' to each other.
#'
#' Warning: for 'negative weight' phrases this method does not do well since it ignores
#' negative documents.
#'
#' @param result  An textreg.result object or a matrix from make.appearance.matrix
#' @family Phrase Vizualization
#' @export
make.similarity.matrix = function( result ) {

	if ( is.textreg.result( result ) ) {
		mat = make.appearance.matrix( result )
	} else {
		mat = result
	}
		
	hits = sqrt( apply( mat, 2, sum ) )
	hits[hits==0] = 1
	
	if ( length(hits) > 1 ) {
		dg = diag( 1/hits )
	} else {
		dg = 1/hits
	}
	sim.mat = dg %*% t(mat) %*% mat %*% dg
	rownames(sim.mat) = colnames(sim.mat) = colnames(mat)
	dim( sim.mat )
	#diag( sim.mat ) = NA
	summary( as.numeric( sim.mat ) )
	stopifnot( max( sim.mat, na.rm=TRUE ) <= 1.000001 )
	sim.mat[ sim.mat > 1 ] = 1
	sim.mat
}



#' Cluster phrases based on similarity of appearance.
#'
#' @description
#' Cluster phrases based on similarity of their appearance in the positive documents.
#' Can also plot this if so desired.
#'
#' Uses hclust() with the ``ward.D'' method on 1-S with S from make.similarity.matrix
#'
#' Warning: for 'negative weight' phrases this method does not do well since it ignores
#' negative documents.
#'
#' @param result A similarity matrix from make.similarity.matrix call or an textreg.result object
#' @param num.groups Number of groups to box.
#' @param plot  Actually plot clustering or just calculate it.
#' @param yaxt  Whether to include a y-axis
#' @param ylab  Label for y-axis
#' @param sub   Subtitle for plot
#' @param main  Title of plot.
#' @param ...   Extra arguments to pass to the plot command.  See par.
#' @family Phrase Vizualization
#' @export
cluster.phrases = function( result, num.groups=5, plot=TRUE, yaxt="n", ylab="", sub="", main="Association of Phrases", ... ) {
	if ( is.textreg.result( result ) ) {
		mat = make.appearance.matrix ( result )
		sim.mat = make.similarity.matrix( result )
	} else {
		mat = NULL
		sim.mat = result
	}
	
	d = 1 - sim.mat[-1,-1]
	dim( d )
	d = as.dist( d )
	d
	
	# or
	# d = dist( t(mat), method="euclidean" )
	
	fit <- hclust(d, method="ward.D")
	
	if( num.groups > nrow( sim.mat ) - 1 ) {
		stop( gettextf( "You cannot form %d groups when you only have %d features", num.groups, length(result$rules) ) )
	}
	
	groups <- cutree(fit, k=num.groups) # cut tree into k clusters
	groups

	# see clustering
	if ( plot ) {
		par( mfrow=c(1,1) )
		tots = apply( mat, 2, sum )
		labs = paste( attr(d,"Labels"), " (", tots[-1], ")", sep="" )
		plot( fit, labels=labs, yaxt=yaxt, ylab=ylab,sub=sub,main=main, ... )
	
		# draw dendogram with red borders around the 5 clusters 
		rect.hclust(fit, k=num.groups, border="red")
	}
	
	# reorder sim matrix to correspond to clustering
	colnames(sim.mat[,-1])[ fit$order ]
	sim.mat = sim.mat[ c(1,1+fit$order), c(1,1+fit$order) ]
	if ( !is.null( mat ) ) {
		mat = mat[ , c(1,1+fit$order) ]
		stopifnot( all( rownames( sim.mat ) == colnames(mat) ) )
	}
	groups = c( NA, groups[fit$order] )
	names(groups)[1] = "*intercept*"
	invisible( list( mat=mat, sim.mat=sim.mat, fit=fit, groups=groups) )
}


#' Generate visualization of phrase overlap.
#'
#' Make simple chart showing which phrases have substantial overlap with other phrases.
#' 
#' @param result   textreg.result object or a similarity matrix from a make.similarity.matrix call.
#' @param count Display counts rather than similarity scores.
#' @param num.groups Number of groups to box.
#' @param use.corrplot  Use the corrplot package of Taiyun Wei (will need to install it).
#' @param ... Extra arguments to pass to the image() plotting command.  See par.
#' @family Phrase Vizualization
#' @export
make.phrase.correlation.chart = function( result, count=FALSE, num.groups=5, use.corrplot=FALSE, ... ) {
	
	lst = cluster.phrases( result, num.groups=num.groups, plot=FALSE )
		
	nlevel = 8
	if ( count || use.corrplot ) {
		img.mat = t(lst$mat) %*% lst$mat
		img.mat = img.mat[-1,-1]
		breaks = c( 0, seq( 1-0.000001,max(img.mat,na.rm=TRUE),length.out= nlevel ) )
	} else {
		img.mat = lst$sim.mat[-1,-1]
		breaks = c( 0, seq(min(img.mat[img.mat>0]),1+0.000001,length.out= nlevel ) )
	}
	
	
	gps = which( lst$groups[-1] != lst$groups[-length(lst$groups)] )
	#gps = c( 1,gps )

	if ( use.corrplot ) {
 		if (!requireNamespace("corrplot", quietly = TRUE)) {
			stop( "need to install corrplot to use the use.corrplot=TRUE option" )
		}

		#if ( count ) {
		#	addCoef.col="black"
		#} else {
		#	addCoef.col=NULL
		#}
		
		corrplot::corrplot( img.mat, is.corr=FALSE, type="lower", method="number", col="black",tl.pos="tp", cl.pos="n" )
		corrplot::corrplot( lst$sim.mat[-1,-1], add=TRUE, is.corr=FALSE, type="upper", diag=FALSE, tl.pos="n", ... )
		#addCoef.col=addCoef.col, ... )
	#corrplot(M,add=TRUE, type="lower", method="number",order="AOE", col="black",
	#diag=FALSE,tl.pos="n", cl.pos="n")
	
		abline( h=(length(lst$groups)-gps) + 0.5, col="red", lty=3 )
		abline( v=gps -0.5, col="red", lty=3 )
		
	} else {
		np = nrow(img.mat)
		par( mar=c(10,10,0.3,0.3), mgp=c(2,1,0) )
		image( 1:np, 1:np, img.mat, xaxt="n", yaxt="n", xlab="", ylab="", 
					breaks=breaks, col = grey(seq(1,0,length.out= nlevel)), ... )
		axis( 1, at=1:nrow(img.mat), labels=rownames(img.mat), cex.axis=0.8, las=2 )
		axis( 2, at=1:nrow(img.mat), labels=rownames(img.mat), cex.axis=0.8, las=2 )
	
		if ( count ) {
			nonz = which( img.mat != 0 )
			corx = 1 + (nonz-1) %% ncol(img.mat) 
			cory = 1 + (nonz-1) %/% ncol(img.mat)
			
			text( corx, cory, as.character( img.mat[nonz] ), col="red", ... )
		}
		abline( h=gps -0.5, col="red", lty=3 )
		abline( v=gps -0.5, col="red", lty=3 )
	}
	
	lst$img.mat = img.mat
	
	invisible( lst )
}




if ( FALSE ) {
	
	## 
	## Demonstrating use of the above
	##
	
	
	##
	## generate and plot clustering
	##
	load( "recycling_results" )
	nres = res_no_bale_C3
	
	res = cluster.phrases( nres, num.groups=3 )
	res$groups
	
	##
	## Now plot "correlation structure"
	##
	make.phrase.correlation.chart( nres )
	make.phrase.correlation.chart( nres, count=TRUE )
	
	
	##
	## Now count shared support: do many phrases overlap?
	##
	
	mat = make.appearance.matrix( nres )
	
	sim.mat = make.similarity.matrix( mat )
	dim( sim.mat )
	
	head( mat )
	p1 = mat[, 3 ]
	
	# count how many documents are shared for each phrase pair
	mt = t(mat) %*% mat
	
	# how many pairwise comparisons are there?
	sum( lower.tri( mt )  )
	
	table( mt[ lower.tri( mt ) ]  )

}
