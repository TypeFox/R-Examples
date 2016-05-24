#' Read NeXML files into various R formats
#' 
#' @param x Path to the file to be read in. Or an \code{\link[XML]{XMLInternalDocument-class}} 
#' or \code{\link[XML]{XMLInternalNode-class}} 
#' @param ... Further arguments passed on to \code{\link[XML]{xmlParse}}
#' @import XML
#' @import httr
#' @aliases nexml_read read.nexml 
#' @export nexml_read read.nexml 
#' @examples
#' # file
#' f <- system.file("examples", "trees.xml", package="RNeXML")
#' nexml_read(f)
#' \dontrun{ # may take > 5 s
#' # url
#' url <- "https://raw.githubusercontent.com/ropensci/RNeXML/master/inst/examples/trees.xml"
#' nexml_read(url)
#' # character string of XML
#' str <- paste0(readLines(f), collapse = "")
#' nexml_read(str)
#' # XMLInternalDocument
#' library("httr")
#' library("XML")
#' x <- xmlParse(content(GET(url)))
#' nexml_read(x)
#' # XMLInternalNode
#' nexml_read(xmlRoot(x))
#' }
nexml_read <- function(x, ...) {
  UseMethod("nexml_read")    
}

#' @export
#' @rdname nexml_read
nexml_read.character <- function(x, ...) {
  if (!any(grepl("^https?://", x), 
                XML::isXMLString(x), 
                file.exists(x))) {
    stop("character input must be a URL, xml string or file path", call. = FALSE)
  }
  # handle remote paths using httr::GET
  if (grepl("^https?://", x)) {
    tmp <- GET(x)
    stop_for_status(tmp)
    x <- content(tmp)
  }
  doc <- xmlParse(x, ...)
  output <- as(xmlRoot(doc), "nexml")
  free(doc) # explicitly free the pointers after conversion into S4
  return(output)
}

#' @export
#' @rdname nexml_read
nexml_read.XMLInternalDocument <- function(x, ...) {
  as(xmlRoot(x), "nexml")
}

#' @export
#' @rdname nexml_read
nexml_read.XMLInternalNode <- function(x, ...) {
  as(x, "nexml")
}


setAs("XMLInternalNode", "phylo", function(from)
  as(as(from, "nexml"), "phylo")
)

read.nexml <- nexml_read

