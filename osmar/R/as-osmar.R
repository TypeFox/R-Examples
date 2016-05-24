#' @include as-osmar-elements.R
{}



subclass <- function(obj, subclass) {
  structure(obj, class = c(subclass, class(obj)))
}



### OSM parsing: #####################################################

element_parse <- function(obj, element) {
  elem <- list()
  elem$elements <- getNodeSet(obj, path = sprintf("//%s", element))
  elem$IDs <- sapply(elem$elements, xmlGetAttr, "id")

  subclass(elem, sprintf("%s_parsed", element))
}



osm_parse <- function(x) {
  stopifnot(class(x) %in% c("XMLInternalElementNode",
                            "XMLInternalNode",
                            "XMLAbstractNode" ))

  osm <- list()
  osm$nodes <- element_parse(x, "node")
  osm$ways <- element_parse(x, "way")
  osm$relations <- element_parse(x, "relation")

  subclass(osm, "osm_parsed")
}



### OSMAR object construction: #######################################

NODES_CLASS <- "nodes"
WAYS_CLASS <- "ways"
RELATIONS_CLASS <- "relations"
OSMAR_CLASS <- "osmar"



osmar_elemclass <- function(obj, subclass) {
  stopifnot(all(sapply(obj, class) == "data.frame"))
  subclass(obj, c(subclass, "osmar_element"))
}



osmar_class <- function(obj) {
  stopifnot(length(obj) == 3)
  #stopifnot(sapply(obj,
  # function(k) class(k)[1])==c("NODE", "WAY", "RELATION"))

  subclass(obj, OSMAR_CLASS)
}


is_osmar <- function(obj) {
  OSMAR_CLASS %in% class(obj)
}

are_osmar <- function(objs) {
  all(OSMAR_CLASS %in% sapply(objs, class))
}



#' Convert OSM-XML to an osmar object
#'
#' Convert a given OSM-XML object (as parsed by
#' \code{\link[XML]{xmlParse}}) to an osmar object.
#'
#' @param xml An OSM-XML object
#'
#' @return
#'   A list (with class attribute \code{osmar}) with three elements:
#'
#'   \describe{
#'
#'     \item{\code{nodes}}{A list with two data frames containing the
#'       attributes and tags of the nodes.}
#'
#'     \item{\code{ways}}{A list with three data frames containing the
#'       attributes, tags, and references of the ways.}
#'
#'     \item{\code{relations}}{A list with three data frames
#'       containing the attributes, tags, and references of the
#'       relations.}
#'
#'   }
#'
#' @aliases osmar
#'
#' @examples
#'   file <- system.file("extdata", "kaufstr.xml", package = "osmar")
#'   raw <- readLines(file)
#'   kaufstr <- as_osmar(xmlParse(raw))
#'
#' @export
as_osmar <- function(xml) {
  #stopifnot(osm_check(xml))

  osm_parsed <- osm_parse(xmlRoot(xml))

  osm_data <- extract_data(osm_parsed)
  osm_attr <- extract_attr(osm_parsed)
  osm_ref <- extract_ref(osm_parsed)

  osmar <- list()

  osmar$nodes <- osmar_elemclass(list(attrs = osm_attr$nodeattr,
                                      tags = osm_data$nodedata),
                                 NODES_CLASS)

  osmar$ways <- osmar_elemclass(list(attrs = osm_attr$wayattr,
                                     tags = osm_data$waydata,
                                     refs = osm_ref$wayref),
                                WAYS_CLASS)

  osmar$relations <- osmar_elemclass(list(attrs = osm_attr$relationattr,
                                          tags = osm_data$relationdata,
                                          refs = osm_ref$relationref),
                                     RELATIONS_CLASS)

  osmar_class(osmar)
}

