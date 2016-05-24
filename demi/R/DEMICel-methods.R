#==============================================================================#
# DEMICel-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize.DEMICel
#==============================================================================#

#' Initializes a \code{DEMICel} object
#' 
#' Initializes a \code{DEMICel} object.
#' 
#' @param .Object A DEMICel object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a \code{DEMICel} object
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMICel" <-
function( .Object, ... ) 
{
	.Object <- callNextMethod( .Object, ... );
	.Object;
}#initialize.DEMICel

#' @import methods
setMethod( "initialize", "DEMICel", initialize.DEMICel );
