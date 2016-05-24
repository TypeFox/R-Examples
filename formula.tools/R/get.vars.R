# -----------------------------------------------------------------------------
# METHOD: get.vars
#
#   Retrieves the variable names from various types of R objects such as calls, 
#   expressions, names and formulas.
#   
#   This method is similar to all.vars except it will expand '.' and other
#   special characters in model formula
#   
#   Returns the variables in order of appearance
# -----------------------------------------------------------------------------


#' Get variable (names) from various R objects
#' 
#' \code{get.vars} extracts variable names from various R objects such as
#' formulas, expressions, calls, symbols, etc.  It is very similar to
#' \code{\link[base]{all.vars}} except that all symbols, etc. are interpolated 
#' to the names of variables.
#' 
#' @param x object to extract vars from.
#' @param data data set/list or environment on which the names are defined
#' @param ... arguments passed to subsequent functions
#' 
#' \code{get.vars} and variant get the variables from objects optionally 
#' interpreting on \code{.} on the data.  This is useful, for example, when you 
#' wish to know what data is used based on a given formula.
#' 
#' Methods/functions beginning with \code{.} are not exported
#' 
#' @return character vector of variables names in order that they appear in
#' \code{x}.
#' 
#' @seealso 
#'   \code{\link[base]{all.vars}}
#' 
#' @examples 
#'   get.vars( Species ~ ., iris )
#'   get.vars( quote( Sepal.Length * Sepal.Width ), iris )
#'   
#' @rdname get.vars
#' @name get.vars
#' @export get.vars

setGeneric( 
  'get.vars', function(x, data=NULL, ...) standardGeneric( 'get.vars' ) 
)


# ---------------------------------------------------------------------
# SIGNATURE: formula
#   For this to work correctly we need to treat the lhs and rhs distinctly
#   and merge the results.
#
#   Some edge cases may not work.
#  
# ---------------------------------------------------------------------

#' @rdname get.vars 
#' @aliases get.vars,formula,ANY-method
#' @export
setMethod( 'get.vars', c( 'formula', 'ANY' ) ,
  # get.vars.form <- 
  function(x, data=NULL, ... ) {
    
    vars.lhs <- get.vars( lhs(x), data=data, ... )

    term.rhs <- terms.formula( x, data=data, ... )
    labels   <- attr( term.rhs, 'term.labels' )
    order    <- attr( term.rhs, 'order' )
    vars.rhs <- labels[ order == 1 ]

    unique( c(vars.lhs, vars.rhs)  )
    
  }
)



# ---------------------------------------------------------------------
#' @rdname get.vars
#' @aliases get.vars,call,ANY-method

setMethod( 'get.vars', c( 'call', 'ANY' ), 
  #  get.vars.call <- function(x,data,...) {
  function( x, data=NULL, ... ) {

    term <- terms( x, data=data, ... )
    return(term)
  
  }
)


# ---------------------------------------------------------------------
#' @rdname get.vars
#' @aliases get.vars,expression,missing-method

setMethod( 'get.vars', c( 'expression', 'missing' ) ,
  function( x, data, ... ) all.vars( x, ... ) 
)


# ---------------------------------------------------------------------
#' @rdname get.vars
#' @aliases get.vars,name,ANY-method

setMethod( 'get.vars', c( 'name', 'ANY' ) ,
  function( x, data, ... ) as.character(x) 
)



# ---------------------------------------------------------------------
#' @rdname get.vars
#' @aliases get.vars,ANY,ANY-methods
setMethod( 'get.vars', c( 'ANY', 'ANY' ), 
  function( x, data, ... ) character(0)
)





# ---------------------------------------------------------------------
#' @rdname get.vars
#' @aliases get.vars,NULL,ANY-methods
setMethod( 'get.vars', c( 'NULL', 'ANY' ), 
  function( x, data, ... ) NULL 
)



