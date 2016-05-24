#' @include osm-descriptors.R
{}



osmsource <- function(obj, subclass) {
  structure(obj, class = c(subclass, "osmsource"))
}



get_osm_data <- function(source, ...) {
  UseMethod("get_osm_data")
}



setOldClass("api")
setOldClass("osmosis")
setOldClass("bbox")
setOldClass("element")
setOldClass(c("node", "element"))
setOldClass(c("way", "element"))
setOldClass(c("relation", "element"))



setGeneric("osm_request",
function(source, what, ...) {
  standardGeneric("osm_request")
})

