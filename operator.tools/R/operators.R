#' Return the _names_ of defined operators.
#' 
#' \code{operators} returns the names of defined operators. Argument
#' \code{types} can be used to select operators of a specified type(s) or
#' GROUPING(s).  See Details for specifics.
#' 
#' @param types A character vector with the types of operators to return. The
#' types may one or more of: 'namespace', 'component', 'indexing', 'sequence',
#' 'arithmetic', 'relational', 'logical', 'tilde', 'assignment', 'help',
#' 'user', or user-defined type specified in a call to
#' \code{\link{setOperator}}.  It may also be one of the special groups:
#' 'REG(ISTERED)', 'UNREG(ISTERED)', 'SPECIAL', 'ALL'. See Details.
#'
#' \code{operators} provides the \bold{names} of defined operators.  These can
#' be either registered operators (using \code{\link{setOperators}}), or
#' unregistered operators definde by the \code{\%any\%} syntax.
#' 
#' By default, only registered operators are returned. This is purely for
#' performance reasons as an exhausting search for \code{\%any\%} functions is
#' expensive.
#' 
#' See \link[base]{Syntax}.for the core R operators
#' 
#' \code{types} may also be one a special operator groupings: \itemize{ \item
#' REG(ISTERED): (Default). Those registered by \code{\link{setOperators}}
#' \item UNREG(ISTERED): Unregisted operators, requires expensive search.
#' \item ALL: All operators, requires expensive search of environments.  \item
#' SPECIAL: All operators defined using the \code{\%any\%} syntax.  }
#' 
#'   
#' @return character vector of unique operator names.
#' 
#' @note The right arrow assignment operators, \code{->} and \code{->>} is not
#' an operator but a syntatic variant.  Consequently, it does not behave
#' properly as an operator.  They are omitted from the operator list as they
#' are not correctly identified as primitives or functions by the R language.
#' 
#' @author Christopher Brown
#' 
#' @seealso \code{\link{Syntax}}, \code{\link{setOperator}},
#' \code{\link{setOperators}}, and the help files on the individial operators.
#' 
#' @references 
#'    \url{http://cran.r-project.org/doc/manuals/R-lang.html}
#'    \url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14310}
#' 
#' @keywords utilities
#' 
#' @examples
#' 
#'  \dontrun{ 
#'   operators()
#'   operators( types="arithmetic" )
#'   operators( types=c("arithmetic","logical" ) )
#'   operators( types='ALL' )
#'   operators( types='REG' )
#'   operators( types='UNREG' )
#'   operators( types='SPECIAL' )
#'  }
#'
#' @export
operators <- function( types="REGISTERED" ) {

  core <- names( .Options$operators )
  
# DEPRECATED: Default was types=NULL, now types=="CORE" 
# CASE: types == NULL --> "ALL" OPERATORS
#  - We test for this most common special case for speed.
  if( is.null(types) ) 
    return( unique( c( core, apropos("^%.*%" ) ) ) )
  
  if( length(types) == 1 && types %in% c("REG","REGISTERED") ) return(core)  

  op.types <- unique( sapply( .Options$operators, function(x) x$type ) )
  ret <- NULL
  for( type in types ) {

    if( type == 'ALL' ) {
      ret <- c(ret, core, apropos("^%.*%" ) ) 
    
    } else if ( type == 'SPECIAL' ) {
      ret <- c(ret, apropos( "^%.*%" ) )

    } else if ( type %in% c( 'UNREG', 'UNREGISTERED' ) ) { 
      special <- apropos( "^%.*%" )
      ret <- c(ret, special[ ! special %in% core ])   
    
    } else if ( type %in% op.types ) {
      ret <- c( ret, core[ sapply( .Options$operators, function(x,y) x$type == y, type ) ] )

    } else {
      stop( "operator type: ", type, " is unknown" ) 
    }

  }

  return(unique(ret)) 

}

