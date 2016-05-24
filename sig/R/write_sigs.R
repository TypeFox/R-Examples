#' Write sigs to file
#' 
#' Writes a list of function signatures to a file.
#' 
#' @param x A list of function signatures.  See details.
#' @param file A file path or connection to write the output to (stdout 
#' by default).
#' @param ... passed to \code{toString.siglist}.
#' @return A character vector of the lines that were written to file is
#' invisibly returned.  Mostly invoked for the side effect of writing 
#' function signatures to a file.
#' @details Where \code{x} is an object of class |code{siglist}, the
#' function essentially calls \code{writeLines(tostring(x))}.
#' If the input  is a single function signature (of class \code{sig}), 
#' then it is coerced into a \code{siglist}.  If the input is an  
#' environment or path to a file, then \code{list_sigs} is called on 
#' the input before writing.
#' @examples
#' #Step by step:
#' #First, list some sigs.
#' utils_sigs <- list_sigs(pkg2env(utils))   
#' #Without a file argument, sigs are just printed to the console.          
#' head(write_sigs(utils_sigs))
#' #Write to a file
#' tmpf <- tempfile("sig", fileext = ".R")
#' write_sigs(utils_sigs, tmpf)    
#' \dontrun{
#' Open the file we've just written
#' shell(tmpf, wait = FALSE)
#' }   
#' #Can also list and write in one line.
#' tmpf2 <- tempfile("sig", fileext = ".R") 
#' write_sigs(pkg2env(grDevices), tmpf2)
#' #Single sigs are coerced to siglists
#' write_sigs(sig(stats::var))           
#' @export
write_sigs <- function(x, file = stdout(), ...)
{
  UseMethod("write_sigs")
}
        
#' @rdname write_sigs
#' @method write_sigs default
#' @export
write_sigs.default <- function(x, file = stdout(), ...)
{
  y <- as.siglist(x, ...)
  write_sigs(y, file, ...)
}
           
#' @rdname write_sigs
#' @method write_sigs siglist
#' @export
write_sigs.siglist <- function(x, file = stdout(), ...)
{
  lines <- toString(x, ...)
  writeLines(lines, file)
  invisible(lines)
}

#write_sigs.sig <- function(x, file, ...)
#{
#  x <- as.siglist(x, ...)
#  NextMethod("write_sigs")
#}
                                   
#' @rdname write_sigs
#' @method write_sigs environment
#' @export
write_sigs.environment <- function(x, file = stdout(), ...)
{
  sig_list <- list_sigs(x, ...)
  write_sigs(sig_list, file = file, ...)
}
                           
#' @rdname write_sigs
#' @method write_sigs character
#' @export
write_sigs.character <- function(x, file = stdout(), ...)
{
  sig_list <- list_sigs(x, ...)
  write_sigs(sig_list, file = file, ...)
}
