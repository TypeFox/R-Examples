# ------------------------------------------------------------------------------
# COERCION FUNCTIONS:
#   fun2name, name2fun, char2fun, fun2char, name2char, char2name
#   as.character.name, as.character.function
#   as.function.character, as.function.name
#   as.name.character, as.name.function
#
#   Functions can be represented as a string, a name in the formal sense or 
#   a definition, i.e. the , function(body) or string   
#   
#   functions 
#   TODO:
#    - Move to function.tools 
#    - recast as as.character, as.name, as.digest, as.s
#
#   See Also: args, body and formals.  
#    - body does not seem to work.  body( `<` ) and body( eval(`<`) ) are both
#      null.
#
# ------------------------------------------------------------------------------

#' Convert between a function and its name and vice versa.
#' 
#' \code{fun2name} compares a function (body) to all defined functions. If an
#' identical match is found to a defined function, that function is returned.
#' NB. This does not search through S4 methods.
#' 
#' \code{name2fun} simply converts its argument to a name and than evals it 
#' to produce a function definition
#' 
#' @param f function 
#' @param x name; more specifically, an object to be converted into a name and eval'd
#' 
#' \code{fun2name} compares the function against existing functions using
#' \code{\link{identical}}. If a match is found, the name of the matching 
#' function ( expressed as a \code{character} ) is returned.
#' 
#' \code{fun2name} will not work for S4 Methods.   
#' 
#' @return 
#'   fun2name: character (name of function)  
#'   name2fun: function
#'  
#' @export
fun2name <- function(f) {
  
  nms <- apropos( ",*", mode="function" )
  
  for( nm in nms ) 
    if( identical( f, eval( as.name(nm) ) ) ) return(nm)

  return(NULL)
    
}

#' @rdname fun2name
#' @export
name2fun <- function(x) eval(as.name(x))

# char2fun(x) name2fun(x)

# as.character.name <- function(x, ... )  HANDLED
# as.character.function <- function(x, ... ) -x is this the same as name.

# as.name.function
# as.name.




