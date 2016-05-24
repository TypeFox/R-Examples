

#' Element descriptors
#'
#' For getting OSM data and finding elements in an \code{\link{osmar}}
#' object one needs to describe the data---here we provide a simple
#' description language.
#'
#' @param object The descriptor; see details
#'
#' @seealso \code{\link{bbox}}
#' @aliases osm_descriptors
#' @rdname osm_descriptors
#' @export
node <- function(object) {
  UseMethod("node")
}

#' @rdname osm_descriptors
#' @export
way <- function(object) {
  UseMethod("way")
}

#' @rdname osm_descriptors
#' @export
relation <- function(object) {
  UseMethod("relation")
}




### Description by ID: ###############################################

elem_by_id <- function(id, subclass) {
  structure(c(id = id), element = subclass, class = c(subclass, "element"))
}


#' @examples
#'   ## Description by ID (*.default):
#'   node(1)
#' @method node default
#' @rdname osm_descriptors
#' @S3method node default
node.default <- function(object) {
  elem_by_id(object, "node")
}

#' @examples
#'   way(1)
#' @method way default
#' @rdname osm_descriptors
#' @S3method way default
way.default <- function(object) {
  elem_by_id(object, "way")
}


#' @examples
#'   relation(1)
#' @method relation default
#' @rdname osm_descriptors
#' @S3method relation default
relation.default <- function(object) {
   elem_by_id(object, "relation")
}



### Description by condition: ########################################


#' @param condition Condition to describe the object
#' @rdname osm_descriptors
#' @export
attrs <- function(condition) {
  structure(list(condition = substitute(condition)), what = "attrs",
            class = "condition")
}

#' @rdname osm_descriptors
#' @export
tags <- function(condition) {
  structure(list(condition = substitute(condition)), what = "tags",
            class = "condition")
}

#' @rdname osm_descriptors
#' @export
refs <- function(condition) {
  structure(list(condition = substitute(condition)), what = "refs",
            class = "condition")
}


#' @examples
#'   ## Description by condition (*.condition):
#'   node(tags(v == "Marienplatz"))
#' @method relation condition
#' @rdname osm_descriptors
#' @S3method node condition
node.condition <- function(object) {
  structure(object, element = "node")
}

#' @examples
#'   ## Description by condition (*.condition):
#'   way(attrs(id == 17458))
#' @method relation condition
#' @rdname osm_descriptors
#' @S3method way condition
way.condition <- function(object) {
  structure(object, element = "way")
}

#' @method relation condition
#' @rdname osm_descriptors
#' @S3method relation condition
relation.condition <- function(object) {
  structure(object, element = "relation")
}

