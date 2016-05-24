#==============================================================================#
# DEMIResult-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize
# getResult
# getGroup
# makeDEMIResultsTable
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMIClust initialization:
#------------------------------------------------------------------------------#

#' Initializes the \code{DEMIResult} object
#' 
#' Initializes the \code{DEMIResult} object.
#' 
#' @param .Object A DEMIResult object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a 'DEMIResult' object.
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMIResult" <-
function( .Object, ... ) 
{
	.Object <- callNextMethod( .Object, ... );
	.Object;
}#initialize.DEMIClust

#' @import methods
setMethod( "initialize", "DEMIResult", initialize.DEMIResult );


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIResult get functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname getResult-methods
#' @aliases getResult,DEMIResult-method
#' @import methods
setMethod( "getResult", signature( object = "DEMIResult" ),
		function( object ) object@result
)#getResult

#' @rdname getGroup-methods
#' @aliases getGroup,DEMIResult-method
#' @import methods
setMethod( "getGroup", signature( object = "DEMIResult" ),
		function( object ) object@group
)#getGroup

#' Returns a \code{data.frame} of the differential expression results
#' 
#' Returns a \code{data.frame} of the differential expression results stored in the \code{DEMIResult}
#' object. It is used internally by DEMI methods.
#' 
#' @param input A \code{list}. Represents a \code{list} of \code{DEMIResult} object.
#' @return Returns a \code{data.frame} of the differential expression analysis results stored in the
#' 		   \code{DEMIResult} objects.
#' @seealso \code{DEMIResult}
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname makeDEMIResultsTable-methods
"makeDEMIResultsTable" <- function( input = "list" )
{
	output <- NULL;
	if ( length( input ) > 0 ) {
		for( i in 1:length( input ) ) {
			if ( class( input[[i]] ) == "DEMIResult" ) {
				demiResult <- getResult( input[[i]] );
				for ( ii in 1:length( demiResult ) ) {
					clusterID <- names( demiResult )[ii];
					output <- rbind( output, ( data.frame( clusterID, demiResult[[clusterID]] ) ) );
				}
			} else {
				#stop( "Error: the results list can only contain objects of class 'DEMIResult'\n" );
				stop( DEMIMessages$DEMIResults$resultsContain );
			}
		}
	} else {
		#stop( "Error: no results have been defined\n" );
		stop( DEMIMessages$DEMIResults$resultsUndefined );
	}
	output <- output[ order( output$FDR ), ];
	return( output );
}#makeDEMIResultTable

