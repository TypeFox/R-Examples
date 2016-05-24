#' @rdname str_uppercase
#' @examples  
#'   str_decapitalize( "ABC" )     # abc
#' @aliases lowercase lower_case
#' @export

  str_lowercase <- function(string) {
      if (!is.atomic(string)) 
      stop("String must be an atomic vector", call. = FALSE)
    if (!is.character(string)) 
      string <- as.character(string)
    
    base::tolower(string)
  }

  # Alternative form:
  # make_str_replace( '([A-Z])', '\\L\\1' ) 

#' @rdname str_uppercase
#' @export
  str_lower_case <- str_lowercase

#' @rdname str_uppercase
#' @export
  str_lower <- str_lowercase

#' @rdname str_uppercase
#' @export
  str_decapitalize <- str_lowercase


#' @rdname str_uppercase 
#' @examples 
#' # is_lower_case 
#'   is_lowercase( 'ABC123' )      # FALSE
#'   is_lowercase( 'abc123' )      # TRUE 
#'   is_lowercase( 'aB'  )         # FALSE 
#'   is_lowercase( '123' )         # TRUE 
#'   
#' @export   
  is_lowercase <- function( string ) 
    ! grepl( "[A-Z]", string, perl=TRUE )

#' @rdname str_uppercase
#' @export
  is_lower_case <- is_lowercase

#' @rdname str_uppercase
#' @export
  is_lower <- is_lowercase
