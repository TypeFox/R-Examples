#' @include Class-PatternOrCharacter.R 
NULL 

#' Extraction operators for Searchable object
#' 
#' Defines  \code{[}, \code{[[}, and \code{$} for Searchable objects
#' 
#' @param x Searchable object
#' @param i character; pattern with potential match modifiers applied,
#' @param j missing; never specified
#' @param drop For matrices and arrays. If TRUE the result is coerced to the 
#'    lowest possible dimension (see the examples). This only works for 
#'    extracting elements, not for the replacement. See \code{\link[base]{drop}}
#'    for further details.
#' @param ... additional arguments. See \code{\link[base]{Extract}}
#' 
#' @param value replacement value for replacement functions 
#'
#' The methods for searching respect the modifiers applied to both \code{x} and
#' \code{i}.
#' 
#' @section \code{[}, \code{[<-} :
#' 
#' \code{[} and \code{[<-} are used for subsetting and replacing 
#' \strong{zero or more} elemenxts of \code{x}. Used with searchable objects, 
#' these operators differ from normal R operations in the following respects:
#' 
#' \itemize{ 
#' 
#' \item The search returns elements of 
#' the target that matches \strong{ANY} of the search patterns. 
#' 
#' \item Unlike the 
#' its normal behavior, \code{\[} does not guarantee the output to have as many
#' elements as elements to \code{pattern}.
#'  
#' \item \code{[} does not return a Searchable object. It is thought that 
#' the return valuable will not be subsequently searched. It is easy to turn 
#' the results into a Searchable object using \code{searchable} however. 
#'
#' \item Unlike for environments and hashes, no constraints exist for ensuring 
#' uniqueness for names in vectors and lists. These structures may contain 
#' multiple elements with the same name. Normal attempts to extract by name 
#' yield only the \strong{first} element that matches the name. Using a 
#' \code{Searchable} patterns match yields all matching elements.
#' 
#' }
#'  
#'  
#' @return 
#'   The values after the extracting methods have been applied:\cr
#'   \code{\[} returns a subset of \code{x}, but which is not Searchable.  \cr
#'   \code{\[\[} and \code{\$} return a sinlge element of \code{x}  \cr
#'   
#'  @seealso
#'    \code{\link{Searchable}}           \cr
#'    \code{\link[base]{Extract}}        \cr
#'    Match mofiers: \code{\link{fixed}}, \code{\link{regex}}, 
#'      \code{\link{coll}} and \code{\link{ignore.case}}
#'    \code{\link{reverse.lookup}}       \cr
#'    
#'    
#' @examples 
#' 
#'   # ATOMIC VECTORS: 
#'     v <- c( a=1, b=2, B=3, c=4, c2=5 )
#'     sv <- searchable(v)
#'       
#'       
#'   # FLEXIBLY FIND ELEMENTS BY NAME 
#'     sv[ regex('c') ]
#'     sv[ fixed('c') ]
#'
#'     sv[ ignore.case('b') ] 
#'                                                                                                                                                                                                                                                                                                                            
#'
#'   # FLEXIBLY REPLACEMENT ELEMENTS BY NAME  
#'     sv[ regex('c.?') ]   <- "3rd"
#'   
#'   
#'   # SET DEFAULT SEARCH FOR TARGET/OBJECT
#'     sv <- searchable(v, case_insensitive = TRUE )         
#'     sv['b']
#'     sv['B']
#'   
#'     sv <- regex(sv)  
#'     sv['c']  
#'
#'     sv <- ignore.case(sv)    
#'     sv['b']                                                                    
#'     sv['c']                  # st  
#'                                        
#'
#'   # USE ON (RECURSIVE) LISTS:
#'     l <- list( a=1, b=2, c=3 )
#'     sl <- searchable(l)                
#'     sl["b"]
#'     sl[ ignore.case("B") ] 
#'     
#'     
#'   # USE WITH MAGRITTR   
#'    \dontrun{
#'     sl[ "B"  %>% ignore.case ]
#'     "b" %>% sl[.]
#'     "B" %>% ignore.case %>% sl[.]
#'    }
#'    
#'      
#' @aliases extract

# ----------------------------------------------------------------
# EXTRACT 
# ----------------------------------------------------------------
  
#' @rdname extract
#' @export     
  setMethod( '[', c(x='Searchable', i='PatternOrCharacter', j='missing'), 
    function(x,i,j,...) {
       
     # ESCAPE HATCH FOR  'std' matching
       if( .is.basic(x,i) ) return( x@.Data[i] ) 
       
       return( x[ .matches(x,i) ] )
     
    }           
  )
  

# ----------------------------------------------------------------
# REPLACE  
# ----------------------------------------------------------------
  
#' @rdname extract
#' @export  
  setReplaceMethod( '[', c(x='Searchable', i='character', j='missing', value='ANY'), 
    function(x,i,value) {
      
      # ESCAPE HATCH
        if( .is.basic(x,i) ) return( `[<-`(x,i,value) )
      
        wh <- which( .matches(x,i) )
      
      x@.Data[wh] <- value
      return(x)
      
    }                
  )
  
  
# --------------------------------------------------------------
# UTILITIES
# --------------------------------------------------------------

# Resolve patterns 
# 
# Resolves the pattern to use from the search 
#
# @param object Searchable object 
# @param pattern object to use as a pattern 
#
# A pattern can be associated with either argument or both; this ensures
# returns the correct pattern to use
# 

.resolve.patterns <- function( object, pattern ) {
   
    # browser()
    
    if( ! pattern  %>% is('Pattern') ) { 
      str <- pattern 
      pattern <- object %>% pattern
      pattern@.Data <- str
    }
    
    return(pattern)
    
}

# Is this a std R search and should use base R comparisons
# search is standard without modifications
.is.basic <- function(object, pattern)
    ! is( pattern, "Pattern" ) &&                      # uses object 
    object@type == 'std'       &&                      # object is std
    ( is.null( object@options$case_insensitive ) ||    # case sensitve 
      ! object@options$case_insensitive 
    ) ||                                            # -OR- 

    is( pattern, "Pattern" )   &&                      # uses pattern
    pattern@type == 'std'      && 
    ( is.null( pattern@options$case_insensitive ) || # case sensitve 
      ! pattern@options$case_insensitive 
    )

# Return logical indication matching elements for name search
#
# Returns indices for use with [, [<- allowing multiple patterns
#
# @return integer 

.matches <- function( object, pattern ) { 
    
  # browser()
  # Determine pattern to use ...
    pattern <- .resolve.patterns( object, pattern )
  
  # APPLY REVERSE LOOKUP BY INVERTING OBJECT
  #  This does not work for recursive, list-like  objects
    # object <- if( .reverse.lookup(pattern) ) invert(object) else object  
        
  return( 
    detect( str=names(object), pattern=pattern ) 
  )
  
}  

