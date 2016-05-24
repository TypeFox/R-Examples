#' Code Chunk Reference
#' 
#' Reference a code chunk.
#' 
#' @param label The code chunk label.
#' @param prefix Logical: return the prefix (chunk type) along with the rank. 
#' @param link ### Logical: Place an HTML link to the anchor?
#' @param type The type of the referenced chunk.
#' @return A text or url reference to the code chunk.
#'
#' @examples
#' figr("test", type="reference")
#' figr("test")
#' figr("test", TRUE, link=FALSE)
#' if(interactive()){
#'   require(knitr)
#'   opts_knit$set(kfigr.prefix=TRUE, kfigr.link=TRUE)
#'   figr("test")
#' }
#'
#' @export
figr <- function(label, prefix, link, type){
  # get defaults if needed
  if(missing(prefix))
    prefix <- opts_knit$get("kfigr.prefix")  
  if(missing(link))
    link <- opts_knit$get("kfigr.link")
  # check type
  if(missing(type)){
    if(!(label %in% names(get("index", envir=anchorenv))))
      stop("First reference must define the type.")
    else
      type <- get("index", envir=anchorenv)[[label]]
  }
  # get reference
  if(prefix)
    txt <- paste(type, index(label, type))
  else
    txt <- index(label, type)
  if(link)
    paste('<a href="#', label,  '">', txt, "</a>", sep="")
  else
    txt
}
