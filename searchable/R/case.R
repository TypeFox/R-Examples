#' Turn on/off case sensitivity for Searchable and Pattern objects 
#' 
#' Functions for affecting the case sensitivity of matching/  
#' 
#' @param object search pattern or target
#' @param ... additional arguments
#' 
#' \code{ignore.case}/\code{case_insensitive} and 
#' \code{use.case}/\code{case_sensitive} control the case sensitivity of the 
#' matching
#' 
#' The default is to perform case-sensitive matching. 
#' 
#' @seealso 
#'   \code{stri_detect_*} from the \code{stringi} package
#'   
#' @examples 
#'   use.case("pattern")     # case-sensitive (Default)
#'   ignore.case("pattern")  # case-insensitive 
#'   
#' @aliases ignore.case use.case case_insensitive case_sensitive

#' @rdname case
#' @export

  ignore.case <- function(object) UseMethod('ignore.case')

#' @rdname case
#' @export
  ignore.case.SearchableOrPattern <- function(object, ...) { 
    object@options$case_insensitive = TRUE 
    return(object)
  }

#' @rdname case
#' @export
  ignore.case.character <- function(object) object %>% std( case_insensitive = TRUE ) 
  
#' @rdname case
#' @export
  ignore.case.default <- function(object) object %>% as.character %>% ignore.case 
  

#' @rdname case
#' @export
  case_insensitive <- function(object) ignore.case(object)


# --------------------------------------------------------------

#' @rdname case
#' @export
  use.case <- function(object) UseMethod('use.case')

#' @rdname case
#' @export
  use.case.SearchableOrPattern <- function(object) { 
    object@options$case_insensitive = FALSE 
    return(object)
  }

#' @rdname case
#' @export
  use.case.character <- function(object) object %>% std( case_insensitive = FALSE )
  
#' @rdname case
#' @export
  use.case.default <- function(object) object %>% as.character %>% use.case 
  

#' @rdname case
#' @export
  case_sensitive <- function(object) use.case(object)
