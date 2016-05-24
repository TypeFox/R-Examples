#' searchable 
#' 
#' Tools For Custom Searches / Subsets / Slices of Named R Objects
#' 
#' The 'searchable' package provides flexibile methods for subseting named 
#' object by matching the names using case (in)sensitivity, regular or 
#' fixed expressions. searches uses the standard '['
#' operator and allows specification of default search behavior to either the
#' search target (named object) and/or the search pattern.
#' 
#' It was designed to make flexible, high performance dictionary and 
#' thesaurus structures.  
#' 
#' @references 
#'   \url{http://stackoverflow.com/questions/5671719/case-insensitive-search-of-a-list-in-r} \cr
#'   \url{http://stackoverflow.com/questions/27085620/which-command-in-r-with-case-insensitive} \cr
#'   \url{http://stackoverflow.com/questions/21450925/r-grep-usage-finding-the-averages} \cr
#'   
#' @seealso
#'   \code{\link{searchable}} \cr
#'   \url{http://cran.r-project.org/web/packages/qdap}
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
#'     sv
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
#'    
#' @docType package
#' @name searchable-package
#' @import magrittr methods stringi

NULL
