#' Use standard matching 
#' 
#' Functions for matching using standard, defaut matching 
#' 
#' @param str search target
#' @param pattern pattern to attempt
#' @param ... supplementary arguments passed to the underlying functions, 
#'        including additional settings for \code{stri_opts_std}
#' @param opts_std list; optional arguments used by stri_*_std functions        
#' @param case_insensitive logical; enable simple case insensitive matching
#'
#' \code{stri_detect_std} is equivalent to \code{str \%in\% pattern} and is 
#' created to provide a parallel to other search methods.
#'  
#' \code{stri_opts_std}
#' 
#' @return 
#'   logical indicating the matching elements in \code{str}
#'  
#' @seealso 
#'   \code{\link[stringi]{stri_detect}}
#' 
#' @examples 
#'   stri_detect_std( letters[1:5], letters[1:2] )  # TRUE TRUE ...
#'   stri_detect_std( letters[1:5], LETTERS[1:2] )  # ALL FALSE 
#'   stri_detect_std( letters[1:5], LETTERS[1:2], opts_std = list(case_insensitive = TRUE ) )
#'   
#' @rdname stri_std            
#' @export

stri_detect_std <- function( str, pattern, ..., opts_std = NULL ) {
 
  if( ! missing(...) )
    opts_std = as.list(opts_std, ... ) 
  
  if( ! is.null(opts_std$case_insensitive) && opts_std$case_insensitive ) 
    toupper(str) %in% toupper(pattern) else 
    str %in% pattern 
  
}
   
  

#' @rdname stri_std            
#' @export
stri_opts_std <- function( case_insensitive = FALSE, ...) { 

  opts <- list(...)
  if ( !missing(case_insensitive) ) opts["case_insensitive"] <- case_insensitive
  opts

}
