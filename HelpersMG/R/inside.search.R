#' inside.search Search a string within files of a folder
#' @title Search a string within files of a folder
#' @author Marc Girondot
#' @return Return an invisible vector with filenames in which the pattern occurs
#' @param path Path of the folder to search in
#' @param pattern Pattern for file names to search in
#' @param text Text to search in files
#' @param showallfilenames logical. Show all the filenames search for in
#' @param fixed logical. If TRUE, pattern is a string to be matched as is. Overrides all conflicting arguments (see gsub)
#' @param ignore.case logical. if FALSE, the pattern matching for text is case sensitive and if TRUE, case is ignored during matching.
#' @param ... Options for readLines(), example warn = FALSE
#' @description Search for a string inside the files of a folder and return where the string is found.\cr
#' The pattern for files that must be included uses regex for filtering.\cr
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' # Search for files in path with names based on pattern that have the string search inside.
#' inside.search(path=".", pattern="*\\.R$", search="embryogrowth")
#' }
#' @export


inside.search <- function(path=".", pattern="*\\.R$", showallfilenames=FALSE, ..., 
                          fixed=TRUE, ignore.case = FALSE, 
                         text=stop("A text to be searched for is necessary")) {
  
  ls <- list.files(path=path, pattern=pattern)
  if (identical(ls, character(0))) {
    warning("No file match this pattern at this path")
    return(invisible(NULL))
  }
  
  filesok <- NULL
  
  for (f in ls) {
 
    fc <- readLines(con = f, ...)
    x <- gsub(pattern=text, replacement="", x=fc, fixed=fixed, ignore.case = ignore.case)
    if (any(x != fc)) {
      cat("file ", f, "\n")
      ll <- c(1:length(fc))[x != fc]
      if (length(ll)>1) {
        cat("lines ", ll, "\n")
      } else {
        cat("line ", ll, "\n")
      }
      filesok <- c(filesok, f)
    } else {
      if (showallfilenames) {
        cat("file ", f, "\n")
        cat("none", "\n")
      }
    }
  }
  return(invisible(filesok))
}
