#' str_uppercase, str_lowercases
#' 
#' convert string to upper or lower case
#' 
#' @param string character; argument to be converted
#' 
#' @section Upper Case and All Caps:
#'   
#'   \strong{uppercase}, \strong{upper_case}, \strong{allcaps}, 
#'   \strong{all_caps} are all synonyms.
#'   
#'   \strong{upper case} or \code{all caps} contains all capital letters. 
#' 
#' @section Lower Case:
#'   \code{str_all_caps} is a synonym for \code{str_uppercase}. 
#'   \code{is_all_caps} is a synonym for \code{is_upper_case}.
#' 
#' @examples
#' # str_uppercase   
#'   str_uppercase( "one flew over the cuckoo's nest" )
#'   str_uppercase( "catch-22" )  
#' 
#'   str_capitalize( "abc" )               # ABC
#'   str_all_caps( "abc" )                 # ABC
#'   
#' @rdname str_uppercase
#' @aliases uppercase upper_case
#' @export

  str_uppercase <- function(string) {
    if (!is.atomic(string)) 
      stop("String must be an atomic vector", call. = FALSE)
    if (!is.character(string)) 
      string <- as.character(string)
    
    base::toupper(string)
  }

  # Alternative form:
  # make_str_replace( '([a-z])', '\\U\\1' ) 

#' @rdname str_uppercase 
#' @export
  str_upper_case <- str_uppercase

#' @rdname str_uppercase 
#' @export
  str_upper <- str_uppercase

#' @rdname str_uppercase 
#' @export
  str_all_caps <- str_uppercase

#' @rdname str_uppercase 
#' @export
  str_allcaps <- str_uppercase

#' @rdname str_uppercase
#' @export
  str_capitalize <- str_uppercase
  

#' @rdname str_uppercase
#' @seealso 
#'   \code{\link{str_is}} 
#' @examples 
#' # is_uppercase 
#'   is_uppercase( 'ABC123' )      # TRUE
#'   is_uppercase( 'abc123' )      # FALSE 
#'   is_uppercase( 'aB'  )         # FALSE 
#'   is_uppercase( '123' )         # TRUE 
#'    
#' @export   
  is_uppercase <- function( string ) {
    if (!is.atomic(string)) 
      stop("String must be an atomic vector", call. = FALSE)
    if (!is.character(string)) 
      string <- as.character(string)
    ! grepl( "[a-z]", string, perl=TRUE )
  }

#' @rdname str_uppercase
#' @export
  is_upper_case <- is_uppercase

#' @rdname str_uppercase
#' @export
  is_upper <- is_uppercase

#' @rdname str_uppercase
#' @export
  is_all_caps <- is_uppercase

#' @rdname str_uppercase
#' @export
  is_allcaps <- is_uppercase
