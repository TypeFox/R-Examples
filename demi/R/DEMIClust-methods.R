#==============================================================================#
# DEMIClust-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize.DEMIClust:
# validDEMIClust:
# getGroup:
# getExperiment:
# getClustMethod:
# getCutoffPvalue:
# getCluster:
# customObject:
# createGroup:
# cluster:
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMIClust initialization:
#------------------------------------------------------------------------------#

#' Initializes the \code{DEMIClust} object
#' 
#' Initializes the \code{DEMIClust} object.
#' 
#' @param .Object A \code{DEMIClust} object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a \code{DEMIClust} object.
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMIClust" <-
function( .Object, ... ) 
{
	.Object <- callNextMethod( .Object, ... );
	#	if the probe clusters need to be created
	if ( length( .Object@cluster ) == 0 ) {
		.Object <- createGroup(.Object);
		.Object <- cluster(.Object);
	}
	.Object;
}#initialize.DEMIClust

#' @import methods
setMethod( "initialize", "DEMIClust", initialize.DEMIClust );


#------------------------------------------------------------------------------#
# DEMIExperiment validation:
#------------------------------------------------------------------------------#

#' Validates the \code{DEMIClust} object
#' 
#' Validates the \code{DEMIClust} object.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns a validated \code{DEMIClust} object.
#' @author Sten Ilmjarv
"validDEMIClust" <-
function( object )
{
	msg <- NULL;
	
	# check that the 'cutoff.pvalue' is numeric
	if ( is.numeric( object@cutoff.pvalue ) == FALSE ) {
		#msg <- paste( msg, "\tError:", sQuote( "cutoff.pvalue" ), "is not numeric\n" );
		msg <- DEMIMessages$parameterNotOfClass( "cutoff.pvalue", "numeric" );
	} else if ( is.numeric( object@cutoff.pvalue ) == TRUE ) {
		if ( object@cutoff.pvalue > 1 || object@cutoff.pvalue < 0 ) {
			#msg <- paste( msg, "\tError:", sQuote( "cutoff.pvalue" ), "has to be set between 0 and 1\n" );
			msg <- DEMIMessages$hasToBeNumericBetween( "cutoff.pvalue", 0, 1 );
		}
	}
	
	# do the verdict
	if ( is.null( msg ) ) TRUE else paste( "\n", msg )
}#validDEMIExperiment

#' @import methods
setValidity( "DEMIClust", validDEMIClust )#setValidity


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIClust get functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname getGroup-methods
#' @aliases getGroup,DEMIClust-method
#' @import methods
setMethod( "getGroup", signature( object = "DEMIClust" ),
		function( object ) object@group
)#getGroup


#' @rdname getExperiment-methods
#' @aliases getExperiment,DEMIClust-method
#' @import methods
setMethod( "getExperiment", signature( object = "DEMIClust" ),
		function( object ) object@experiment
)#getExperiment

#' @rdname getClustMethod-methods
#' @aliases getClustMethod,DEMIClust-method
#' @import methods
setMethod( "getClustMethod", signature( object = "DEMIClust" ),
		function( object ) object@clust.method
)#getClustMethod

#' @rdname getCutoffPvalue-methods
#' @aliases getCutoffPvalue,DEMIClust-method
#' @import methods
setMethod( "getCutoffPvalue", signature( object = "DEMIClust" ),
		function( object ) object@cutoff.pvalue
)#getCutoffPvalue

#' @rdname getCluster-methods
#' @aliases getCluster,DEMIClust-method
#' @import methods
setMethod( "getCluster", signature( object = "DEMIClust" ),
		function( object ) object@cluster
)#getClusters


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIClust other functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname customObject-methods
#' @aliases customObject,DEMIClust-method
#' @import methods
setMethod( "customObject", signature( object = "DEMIClust" ),
		function( object ) {
			#	the groupA and groupB are only set if the 'DEMIClust' object
			#	is created automatically instead of being specified by the user
			if ( length( object@group@groupA > 0 ) & ( length( object@group@groupB ) > 0 ) ) {
				return( FALSE );
			} else {
				return( TRUE );
			}
		}
)#customObject

#' @rdname createGroup-methods
#' @aliases createGroup,DEMIClust-method
#' @import methods
setMethod( "createGroup", signature( object = "DEMIClust" ),
		function( object ) {
			indexA <- grep( as.character( getGroup( object )@groupA ), colnames( getNormMatrix( object@experiment ) ) );
			indexB <- grep( as.character( getGroup( object )@groupB ), colnames( getNormMatrix( object@experiment ) ) );
			if ( length( indexA ) == 0 ) {
				#stop( paste( "\tThere are no CEL files that include the group name", paste( "'", as.character( getGroup( object )@groupA ), "'", sep = "" ), "\n",
				#				"\tThe", sQuote( "groups" ), "needs to be set with common names that are used in CEL file names to represent common condition\n") );
				stop( DEMIMessages$DEMIClust$noCELFilesWithGroupname( as.character( getGroup( object )@groupA ) ) );
			} else if ( length( indexB ) == 0 ) {
				#stop( paste( "\tThere are no CEL files that include the group name", paste( "'", as.character( getGroup( object )@groupB ), "'", sep = "" ), "\n",
				#				"\tThe", sQuote( "groups" ), "needs to be set with common names that are used in CEL file names to represent common condition\n") );
				stop( DEMIMessages$DEMIClust$noCELFilesWithGroupname( as.character( getGroup( object )@groupB ) ) );
			}
			object@group <- DEMIGroup( groupA = getGroup( object )@groupA,
									   groupB = getGroup( object )@groupB,
									   indexA = indexA,
									   indexB = indexB );
			return( object );
		}
)#createGroups

#' @rdname cluster-methods
#' @aliases cluster,DEMIClust-method
#' @import methods
setMethod( "cluster", signature( object = "DEMIClust" ),
		function( object ) {
			# check for more then two samples for both groups
			group <- getGroup( object );
			fun <- getClustMethod( object );
			if ( identical( fun, demi.wilcox.test ) == TRUE ) {
				#cat( "# Clustering with demi built-in functions\n" );
				cat( DEMIMessages$DEMIClust$usingBuiltIn );
				if ( length( getIndexA( group ) ) > 2 && length( getIndexB( group ) ) > 2 ) {
					fun <- getClustMethod( object );
					object@cluster <- fun( object );
				} else {
					#cat( "\tWe are using greater or lower comparisons since at least one of the groups contains no more then two files.\n" );
					cat( DEMIMessages$DEMIClust$usingCreaterOrLower );
					object@cluster <- demi.comp.test( object );
				}
			}
			else {
				#cat( "# Clustering with user provided custom function\n" );
				cat( DEMIMessages$DEMIClust$usingUserProvided );
				fun <- getClustMethod( object );
				object@cluster <- fun( object );
			}
			return( object );
		}
)#cluster
