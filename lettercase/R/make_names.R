#' make_names
#' 
#' Make syntactically valid names out from character vectors replacing 
#' \code{.} (dot) with \code{_} underscore.  Hhis is useful when you wish
#' to use snake_case naming convention or are using SQL.
#' 
#' @param names character; vector to be coerced to syntactically valid names. 
#' This is coerced to character if necessary.
#' 
#' @param unique logical; if TRUE, the resulting elements are unique. This may be 
#' desired for, e.g., column names.
#' 
#' @param leading_ What to replace leading '_' and '.' with. Can be only: 
#' A-Z, a-z, ., or "" (Defaults)
#' 
#' Calls \code{nake.names} and then replaces \code{.} by \code{_}
#' See \code{\link[base]{make.names}} for details.
#' 
#' Multiple consecutive underscores are replaced by a single underscore.
#' 
#' Names the end up with leading underscores are replaced with 
#' \code{ leading_ } which can be a string of any length beginning with 
#' a letter or '.'. The default is to drop leading underscores.   
#' 
#' This function is idempotent -- multiple application of the
#' function do not change the results.
#' 
#' @return a character vector containing  
#'
#' @references 
#'   \url{https://en.wikipedia.org/wiki/Snake_case} \cr
#'   \url{http://titlecase.com}
#' 
#' @author Christopher Brown   
#' 
#' @seealso \code{\link[base]{make.names}}, \code{\link[base]{make.unique}}
#' 
#' @examples
#'   make_names(c("foo and bar", "foo-and-bar"), unique = TRUE)
#'   # "foo_and_bar"   "foo_and_bar_1"
#'   
#'   make_names(c("foo and bar", "foo.and_bar"), unique = FALSE)
#'   # "foo.and.bar"  "foo_and_bar"
#'   
#'   make_names(c("foo and bar", "foo.and_bar"), unique = TRUE)
#'   # "foo_and_bar"   "foo_and_bar_1"
#'   
#'   make_names( c(".foo", "_bar") )  # "foo" "bar"
#' 
#'   make_names( c(".foo", "_bar"), leading="." )  # ".foo" ".bar" 
#' 
#' @export

make_names <- function( names, unique=FALSE, leading_ = '' ) { 

   if( ! substr( leading_, 1, 1 ) %in%  c( '', '.', letters, LETTERS ) )
     stop( "'leading_' must be a lower- or uppercase letter or '.'" )
  
   # USE make.names
     names <- sub( '^[^A-Za-z\\.]+', '.', names )# replace leading non-lead character with .
     names <- make.names( names, allow_ = TRUE ) # make.names, allow underscores
     names <- gsub( '\\.', '_', names )          # replace . -> _
     names <- gsub( "\\_+", "_", names )         # replace multiple leading _ with single _    
   
     if( unique ) names <- make.unique( names, sep="_")  # make.unique
 
   # REPLACE LEADING _ WITH leading_ 
     # leading <- grepl( '^_', names )
     # substr( names[ leading ], 1, 1 ) = leading_  
     names <- gsub( '^_', leading_, names)
   
  
   return(names)
  
}
