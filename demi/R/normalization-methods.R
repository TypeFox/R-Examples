#==============================================================================#
# normalization-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# celMatrixNormalize
# norm.rrank
# norm.rma
# norm.quantile
#==============================================================================#

#' @rdname celMatrixNormalize-methods
#' @aliases celMatrixNormalize,DEMIExperiment,function-method
#' @import methods
setMethod( "celMatrixNormalize", signature( object = "DEMIExperiment", fun = "function" ),
		function( object, fun ) {
			#cat( "*Normalizing expression values" );
			cat( DEMIMessages$normalization$main );
			#	load the original expression values
			orgdata = getCelMatrix( object );
			#	keep only probes that match to specified criteria
			blatTable <- getAlignment( object );
			probesInBlat <- unique( blatTable$probeID )
			orgdata <- orgdata[ ( rownames( orgdata ) %in% probesInBlat ) == TRUE, ]
			#	calculate the normalized expression values
			ndata <- fun( orgdata );
			#	add normalized matrix to 'DEMIExperiment' object
			object@exprsData@normMatrix <- as.matrix( ndata );
			return( object );
		}
)#celMatrixNormalize

#' @rdname norm.rrank-methods
#' @aliases norm.rrank,matrix-method
#' @import methods
setMethod ( "norm.rrank", signature( object = "matrix" ),
		function( object ) {
			#cat( " - using 'relative rank' as the normalization method\n" );
			cat( DEMIMessages$normalization$normrrank )
			result <- apply( object, 2, function( x ) { rank( x, ties.method = "max" ) / length( x ) * 100 } );
			return( result );
		}
)#norm.rrank

#' @rdname norm.rrank-methods
#' @aliases norm.rrank,numeric-method
#' @import methods
setMethod ( "norm.rrank", signature( object = "numeric" ),
		function( object ) {
			#cat( " - using 'relative rank' as the normalization method\n" );
			cat( DEMIMessages$normalization$normrrank )
			result <- rank( object, ties.method = "max" ) / length( object ) * 100;
			return( result );
		}
)#norm.rrank

#' @rdname norm.quantile-methods
#' @aliases norm.quantile,matrix-method
#' @import methods
setMethod( "norm.quantile", signature( object = "matrix" ),
		function( object ) {
			#cat( "- using 'quantile normalization' as the normalization method\n" );
			cat( DEMIMessages$normalization$normquantile );
			# create a sorting index by columns in ascending order of values
			data <- object;
			imatrix <- apply( data, 2, order );
			nr <- nrow( data );
			# calculate the row-wise average of sorted values
			rank_av <- apply( imatrix, 1, function( x ) { mean( data[ x + 0:( length( x ) - 1 ) * nr ] ) } );
			# create normalized expression matrix
			result <- apply( imatrix, 2, function( x ) { z <- numeric( length( x ) ); z[ x ] <- rank_av; as.numeric( z ) } );
			rownames( result ) <- rownames( data );
			return( result );
		}
)#norm.quantile
