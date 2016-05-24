#' Work with variables as categorical or continuous
#' 
#' These functions facilitate working with variables as categorical or continous
#' rather than logical, integer, numeric, factor, character, .. 
#' 
#' @param x object
#' @param names logical; whether to return the names of the variables instead of
#' their index
#' @param ... arguments passed to other functions
#' 
#' @description
#' These functions are used to identify which/if a variable or variables are 
#' categorical or continuos.  \code{is.cat} and \code{is.cont} take single 
#' variable arguments.  \code{which.cat} and \code{which.cont} take a list or 
#' data.frame or any structures whose elements have a class method.
#'   
#' @return
#' \code{is.cat} returns \code{TRUE} if \code{x} is \code{character}, 
#' \code{factor} or \code{logical}.  It is \code{FALSE} otherwise.
#' 
#' \code{is.cont} returns \code{TRUE} if \code{x} is \code{numeric}, 
#' \code{integer}, \code{complex}, \code{Date}, \code{POSIXct}
#' 
#' \code{which.cat} and \code{which.cont} report which variables in an object 
#' are categorical and which are continuous.  By default, interger indices are
#' return.  If \code{names=TRUE}, the names of the variables are returned 
#' instead.
#' 
#' @seealso
#'    \code{\link{which}}, \code{\link{is}} 
#'    
#' @examples
#'  \dontrun{
#'   is.cat(letters)          # TRUE
#'   is.cat(factor(letters))  # TRUE
#'   is.cat(TRUE)             # TRUE 
#'   is.cat(FALSE)            # TRUE
#'   is.cat(1:10)             # FALSE
#'   is.cat(rnorm(10))        # FALSE  
#'   is.cat( now() )          # FALSE 
#'   
#'   is.cont(letters)         # FALSE
#'   is.cont(factor(letters)) # FALSE
#'   is.cont(TRUE)            # FALSE 
#'   is.cont(FALSE)           # FALSE
#'   is.cont(1:10)            # TRUE
#'   is.cont(rnorm(10))      # TRUE  
#'   
#'   which.cat(iris)
#'   which.cat( iris, names=TRUE )
#'   
#'   which.cont( iris )
#'   which.cont( iris, names=TRUE )
#'  }
#' @docType methods
#' @rdname catcont


#' @rdname catcont
#' @aliases is.cat
#' @name is.cat
#' @export is.cat
setGeneric( 'is.cat', function(x, ...) standardGeneric( 'is.cat' ) )

#' @rdname catcont
#' @aliases is.cat,character-method
setMethod( 'is.cat', 'character', function(x) TRUE )

#' @rdname catcont
#' @aliases is.cat,factor-method
setMethod( 'is.cat', 'factor',    function(x) TRUE )

#' @rdname catcont
#' @aliases is.cat,logical-method
setMethod( 'is.cat', 'logical',   function(x) TRUE )

#' @rdname catcont
#' @aliases is.cat,ANY-method
setMethod( 'is.cat', 'ANY', function(x) FALSE ) 


# is.cat.default  <- function(x, classes = c( 'character', 'factor', 'logical' ) ) 
#   class( x ) %in% classes       


#' @rdname catcont
#' @aliases is.cont
#' @export  is.cont
#' @name is.cont
setGeneric( 'is.cont', function(x, ...) standardGeneric( 'is.cont' ) )

#' @rdname catcont
#' @aliases is.cont,numeric-method
setMethod( 'is.cont', 'numeric', function(x) TRUE ) 

#' @rdname catcont
#' @aliases is.cont,integer-method
setMethod( 'is.cont', 'integer', function(x) TRUE )

#' @rdname catcont
#' @aliases is.cont,complex-method
setMethod( 'is.cont', 'complex', function(x) TRUE ) 

#' @rdname catcont
#' @aliases is.cont,Date-method
setMethod( 'is.cont', 'Date'   , function(x) TRUE ) 

#' @rdname catcont
#' @aliases is.cont,POSIXct-method
setMethod( 'is.cont', 'POSIXct', function(x) TRUE )


#' @rdname catcont
#' @aliases is.cont,factor-method
setMethod( 'is.cont', 'factor',  function(x) FALSE )  # REQUIRED

#' @rdname catcont
#' @aliases is.cont,ANY-method
setMethod( 'is.cont', 'ANY', function(x) FALSE ) 

# is.cont.default <- function(x, classes = c( 'numeric', 'integer', 'Date' ) )
#   class( x ) %in% classes 



## WHICH ## 


#' @aliases which.cat
#' @rdname catcont
#' @export
which.cat  <- function(x, ..., names = FALSE ) 
{ 
  ret <- which( unlist( lapply( x, is.cat, ... ) ) )
  return( if( names ) names(ret) else ret  )
}


    
#' @aliases which.cont
#' @rdname catcont
#' @export 
which.cont  <- function(x, ..., names = FALSE ) 
{ 
  ret <- which( unlist( lapply( x, is.cont, ... ) ) )
  return( if( names ) names(ret) else ret  )
}
