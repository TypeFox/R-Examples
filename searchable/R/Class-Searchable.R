#' @title Searchable
#' 
#' @description
#' \code{searchable} makes a named object a \code{Searchable} target, optionally
#' specifying the default search options.
#'   
#' @param object searchable object or object to be made searchable
#' @param type character; the type of search to perform 
#' @param ... additional arguments defining the search pattern. See 
#'   \code{?pattern} for details.
#'  
#' @details 
#' 
#' The searchable class allows 'stringr/i'-like searches using \code{\[} and 
#' \code{\[<-} operators. The following search types are supported: 
#' 
#' \itemize{ 
#'   \item \code{std} standard R matching, the default
#'   \item \code{regex} for regular expression matching, 
#'   \item \code{fixed} for fixed string matching, 
#'   \item \code{coll} for collation matching, 
#' }
#' 
#' Class \code{Searchable} objects allow customizations of how R's  
#' \code{\[} operator match objects' names. 
#' 
#' @section Diffences from stringr:
#' 
#' \code{stringr}  and \code{stringi} are general purpose string manipulations 
#' library allowing flexible search and pattern matching against character 
#' strings. The \code{searchable} package applies this type of matching  
#' to objects' names using the standard \code{\[} accessor.  Thus, 
#' 
#'   \code{ searchable(sv)[ regex('b') ]} 
#' 
#' returns objects the subset of whose names contain 'b'.
#'  
#' Unlike \code{stringr/i}, \code{searchable} allows search specification 
#' to applied to either the search pattern or search target. 
#' When applied to the target, a default search method is configured. All 
#' subsequent searches of the searchable target will use this default pattern. 
#'   
#' The search method can be specified with the \code{type} argument of the 
#' \code{searchable} function or any of match-modifying functions, 
#' e.g. \code{fixed}, \code{regex}, \code{coll}, \code{ignore.case}, etc. 
#' See examples. 
#' 
#' When modifiers are applied to both target and pattern, \strong{modifers 
#' applied to the pattern take precedence} and the target's modifiers are 
#' disabled. 
#' 
#' 
#' @section Differences from base R:
#' 
#' \code{searchable} is designed to be minimally invase. When no search types 
#' or options are specified, mathcing defaults to R's normal behavior. 
#' 
#' Here are the other differnece from standard R operations: 
#' 
#' \itemize{
#'  
#'   \item \code{\$} and \item{\[\[} are unaltered by the package. It is 
#'   unclear, how these operators might accommodate the indeterminate number of 
#'   matches.  
#'         
#'   \item Searches using multiple patterns recylce the patterns, but rather 
#'   return elements that match any of the patterns.  
#'   
#'   \item In base R, there is output value every element of input argument, 
#'         \code{i}. Input elements that do not match a named element of 
#'         \code{x} return \code{NA}. Because of the indeterminant number of 
#'         matches given a pattern search against a \code{searchable} object, 
#'         there is no guarantee that a search pattern have a match. If no 
#'         matches are found, a zero-length object is returned. (This may change
#'         to \code{NA} to be more consisitent.) 
#'         
#'   \item Results do not yield a Searchable object, but the superclass that
#'         the searchable class wraps. See \strong{Value} below.
#'          
#' }  
#' 
#' 
#' @section replacement:
#' 
#' \code{searchable} can be used to replace objects as well. See \code{?extract} 
#' for additional exemples.
#' 
#' @section multiple dimension objects:
#' 
#' Multiple dimension ojects such as data.frames, data.tables, matrices and 
#' arrays are not supported at this time.
#' 
#' 
#' @return 
#'   By default, extraction from a searchable objects does not produce a subset 
#'   that is also searchable. It is assumed that in most cases, the developer 
#'   will not want another searchable object and only wish to have the subclass.  
#'   
#' @note 
#'   - Environments cannot be (easily) be made "searchable" due to the way the 
#'     they are implemented.
#'     
#'   - The extraction methods for searchable objects are (at present) limited to
#'     only one pattern. This may change in the future.  
#'   
#' @seealso
#'   \code{\link{extract}}              \cr 
#'   \code{\link[stringi]{stri_detect_regex}}       \cr
#'   \code{\link{reverse.lookup}}       \cr
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
#' @rdname searchable
#' @exportClass Searchable
#' @export
#' @include Class-Searchables.R

   Searchable <- setClass( 'Searchable'
     , contains = 'Searchables' 
     , representation = representation( 'Searchables', type='character', options='list')
     , prototype( NA_character_, type = 'std', options=list() ) 
   )

   # setClass( 'Searchable', contains = 'SearchableOrPattern' )
    

# CONSTRUCTOR
# NB. compare with pattern 
#' @rdname searchable
#' @export
 
  searchable <- function( object, type='std', ... ) { 
    
    # TRAP NON-NAMED OBJECTS
    if( object  %>% attr('names')  %>% is.null ) 
      stop( 'Only objects with a names attribute can be made searchable.')
    
    return( 
      new( 'Searchable', object, type=type, options=list(...) )  
    )
  }  

  

# METHOD: show
#' @rdname searchable

  setMethod('show', 'Searchable', 
    function(object) {
      
      cat( 'searchable ' ) 
      cat( class(object@.Data) )
      cat( ' object using', .describe_pattern(object) )       
      show( object@.Data[ 1:length(object@.Data) ] )  # REMOVE attributes
      invisible(NULL)
      
    }
  )
  