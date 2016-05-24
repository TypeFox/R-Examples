#' Generate matrix describing gradient descent path of textreg.
#'
#' Generate a matrix of the sequence of features as they are introduced with the textreg gradient descent 
#' program along with their coefficients with each step of the descent.
#'
#' @family plot.path.matrix
#' @export
#' @param res A textreg.result object.
#' @examples
#' data( testCorpora )
#' testI = testCorpora$testI
#' res = textreg( testI$corpus, testI$labelI, c("frog","goat","bat"), C=2, verbosity=0 )	
#' make.path.matrix( res )
make.path.matrix = function( res ) {

	stps = res$path[[1]]
	stps.name = res$path[[2]]
	
	fets = unique( stps.name )
	
	mat = matrix( 0, ncol= length(fets), nrow=1+(length(stps)/2) )
	cur = rep( 0, length(fets) )
	names(cur) = fets

	for ( i in 1:(length(stps)/2) ) {
		cur[ stps.name[[i]] ] = cur[ stps.name[[i]] ] + stps[i]
		cur[ stps.name[[i+1]] ] = cur[ stps.name[[i+1]] ] + stps[i+1]
		mat[ i+1, ] = cur
	}
	attributes( mat )$features = fets
	mat
}


#' Plot optimization path of textreg.
#' 
#' Plot the sequence of features as they are introduced with the textreg gradient descent 
#' program.
#'
#' @export
#' @param path.matrix Either a textreg.result object or a matrix from the make.path.matrix call.
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param bty Box for plot
#' @param ... Arguments to be passed to the matplot() command.
#' @family plot.path.matrix
path.matrix.chart = function( path.matrix, xlab="step", ylab="beta", bty="n", ... ) {
	if ( is.textreg.result( path.matrix ) ) {
		path.matrix = make.path.matrix( path.matrix )
	}
	fets = attributes(path.matrix)$features
	
	matplot( (1:nrow(path.matrix))/2, path.matrix, lty=1, col=1:length(fets), xlab=xlab, ylab=ylab, bty=bty, type="n", ...)
	abline( h=0, lty=3, col="grey" )
	matplot( (1:nrow(path.matrix))/2, path.matrix, type="l", lwd=2, lty=1, col=1:length(fets), bty="n", pch=19, add=TRUE )
	
	legend( "topright", fets, lty=1,col=1:length(fets), bty="n", cex=0.8, bg="white")
}



#' Plot the sequence of features as they are introduced with the textreg gradient descent 
#' program.
#' 
#' Simply calls path.matrix.chart.
#' @seealso path.matrix.chart
#'
#' @export
#' @param x A textreg.result object.
#' @param ... Parameters to be passed to path.matrix.chart.
#' @family plot.path.matrix
plot.textreg.result = function( x, ... ) {
	path.matrix.chart( x, ... )
}
