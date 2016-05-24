##newEmptyWiki.r
##2013-07-01 dmontaner@cipf.es
##TiddlyWikiR library

##' @name newEmptyWiki
## @docType function
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords wiki template
##' @seealso \code{file.copy}
##' 
##' @title Creates a TiddlyWiki template to start working.
##' 
##' @description An empty TiddlyWiki file is created.
##' It is intended to be used as a template where you will "write" your report
##' and where you will automatically insert chunks of information
##' taken from your R session results.
##' 
##' @details The file is copied form the examples folder in the TiddlyWikiR library.
##' It may not be the latest version of TiddlyWiki so you may want to visit
##' \url{http://tiddlywiki.com} to get the newest one.
##'
##' The TiddlySaver.jar file is also copied (unless indicated).
##' This file is needed by some browsers to be able to save changes in TiddlyWiki 
##' (see \url{http://tiddlywiki.com/#TiddlySaver}).
##' 
##' Some useful plugins are installed within this local version 
##' (see \url{http://tiddlywiki.com/#Plugins}).
##' 
##' @param file the name of the output file.
##' @param overwrite logical. If TRUE the destination file will be overwritten if it exists.
##' @param TiddlySaver logical. If TRUE the "saver" file will be copied to the same directory as the TiddlyWiki file.
##' 
##' @return The function tries to create an empty TiddlyWiki template file.
##' A value of TRUE or FALSE is returned if the file could be created or not.
##'
##' @examples
##' \dontrun{
##' newEmptyWiki ("myTemplate.html")
##' browseURL ("myTemplate.html")
##' }

##' @export
newEmptyWiki <- function (file, overwrite = FALSE, TiddlySaver = TRUE) {
  sampl.file <- file.path (system.file (package = "TiddlyWikiR"), "examples", "empty.html")
  saver.file <- file.path (system.file (package = "TiddlyWikiR"), "examples", "TiddlySaver.jar")

  out <- file.copy (from = sampl.file, to = file, overwrite = overwrite)
  
  if (TiddlySaver) {
    file.copy (from = saver.file, to = file.path (dirname (file), "TiddlySaver.jar"), overwrite = overwrite)
  }

  return (out)
}
