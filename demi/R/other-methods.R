#==============================================================================#
# other-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cleanorganismname
#==============================================================================#

#' Cleans the organism name from redundant characters
#' 
#' Cleans the organism name from redundant. It is used internally in DEMI
#' analysis.
#' 
#' @param organism A \code{character}. A \code{character} specifing the organism name.
#' @return A \code{character} representing clean organism name.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname cleanorganismname-methods
cleanorganismname <- function( organism )
{
	organism <- gsub( "_", "", organism );
	organism <- gsub( "-", "", organism );
	organism <- gsub( " ", "", organism );
	return( organism );
}
