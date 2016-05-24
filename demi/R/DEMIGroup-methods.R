#==============================================================================#
# DEMIGroup-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize
# getIndexA
# getIndexB
# getGroupA
# getGroupB
# getGroupNames
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMIGroup initialization:
#------------------------------------------------------------------------------#

#' Initializes the \code{DEMIGroup} object
#' 
#' Initializes the \code{DEMIGroup} object.
#' 
#' @param .Object A DEMIGroup object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a 'DEMIGroup' object.
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMIGroup" <-
		function( .Object, ... ) 
{
	.Object <- callNextMethod( .Object, ... );
	.Object;
}#initialize.DEMIGroup

#' @import methods
setMethod( "initialize", "DEMIGroup", initialize.DEMIGroup );


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIGroup get functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname getIndexA-methods
#' @aliases getIndexA,DEMIGroup-method
#' @import methods
setMethod( "getIndexA", signature( object = "DEMIGroup" ),
		function( object ) object@indexA
)#getIndexA

#' @rdname getIndexB-methods
#' @aliases getIndexB,DEMIGroup-method
#' @import methods
setMethod( "getIndexB", signature( object = "DEMIGroup" ),
		function( object ) object@indexB
)#getindexB

#' @rdname getGroupA-methods
#' @aliases getGroupA,DEMIGroup-method
#' @import methods
setMethod( "getGroupA", signature( object = "DEMIGroup" ),
		function( object ) object@groupA
)#getgroupA

#' @rdname getGroupB-methods
#' @aliases getGroupB,DEMIGroup-method
#' @import methods
setMethod( "getGroupB", signature( object = "DEMIGroup" ),
		function( object ) object@groupB
)#getgroupB

#' @rdname getGroupNames-methods
#' @aliases getGroupNames,DEMIGroup-method
#' @import methods
setMethod( "getGroupNames", signature( object = "DEMIGroup" ),
		function( object ) object@groupNames
)#getGroupNames
