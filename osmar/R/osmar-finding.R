#' @include osm-descriptors.R
#' @include osmar-subsetting.R
{}



#' Find element for a given condition
#'
#' @details
#'   The basis of an \code{\link{osmar}} object are
#'   \code{data.frame}s; therefore the \code{condition} principally
#'   follows the rules for \code{\link[base]{subset}}: logical
#'   expression indicating elements or rows to keep.
#'
#'   Furthermore, one has to define on which element and which data
#'   of the \code{\link{osmar}} object the condition applies:
#'   \code{element(data(condition))}, see \code{\link{osm_descriptors}}.
#'
#' @examples
#'   data("muc", package = "osmar")
#'   find(muc, node(tags(v == "Marienplatz")))
#'   find(muc, node(tags(v %agrep% "marienplatz")))
#'   find(muc, node(attrs(id == 19475890)))
#'   find(muc, way(tags(k == "highway" & v == "pedestrian")))
#'
#' @param object An \code{\link{osmar}} object
#' @param condition A condition for the element to find; see details
#'   section.
#'
#' @return The ID of the the element
#'
#' @family finding
#' @seealso binary_grep
#'
#' @export
find <- function(object, condition) {
  stopifnot(is_osmar(object))
  stopifnot(attr(condition, "element") %in% c("node", "way", "relation"))

  handler <- sprintf("find_%s", attr(condition, "element"))
  do.call(handler, list(object, condition))
}



find_node <- function(object, ...) {
  UseMethod("find_node")
}

find_node.osmar <- function(object, ...) {
  find_node.nodes(object$nodes, ...)
}

find_node.nodes <- function(object, condition) {
  #stopifnot(class(condition) == "call")
  stopifnot(attr(condition, "what") %in% c("attrs", "tags"))

  what <- attr(condition, "what")

  ## id <- subset(object[[what]], eval(condition$condition), select = id)$id
  ## ... the above doesn't work when calling find_node.nodes
  ## inside a function; therefore (from subset.data.frame):
  r <- eval(condition$condition, object[[what]], parent.frame(3))
  if (!is.logical(r))
    stop("'subset' must evaluate to logical")
  r <- r & !is.na(r)
  id <- object[[what]]$id[r]

  if ( length(id) == 0 )
    id <- as.numeric(NA)

  id
}



find_way <- function(object, ...) {
  UseMethod("find_way")
}

find_way.osmar <- function(object, ...) {
  find_way.ways(object$ways, ...)
}

find_way.ways <- function(object, condition) {
  stopifnot(attr(condition, "what") %in% c("attrs", "tags", "refs"))

  what <- attr(condition, "what")

  ## id <- subset(object[[what]], eval(condition$condition), select = id)$id
  ## ... see find_nodes.node for the explanation.
  r <- eval(condition$condition, object[[what]], parent.frame(3))
  if (!is.logical(r))
    stop("'subset' must evaluate to logical")
  r <- r & !is.na(r)
  id <- object[[what]]$id[r]

  if ( length(id) == 0 )
    id <- as.numeric(NA)

  id
}



find_relation <- function(object, ...) {
  UseMethod("find_relation")
}

find_relation.osmar <- function(object, ...) {
  find_relation.relations(object$relations, ...)
}

find_relation.relations <- function(object, condition) {
  stopifnot(attr(condition, "what") %in% c("attrs", "tags", "refs"))

  what <- attr(condition, "what")
  ## id <- subset(object[[what]], eval(condition$condition), select = id)$id
  ## ... see find_nodes.node for the explanation.
  r <- eval(condition$condition, object[[what]], parent.frame(3))
  if (!is.logical(r))
    stop("'subset' must evaluate to logical")
  r <- r & !is.na(r)
  id <- object[[what]]$id[r]

  if ( length(id) == 0 )
    id <- numeric(NA)

  id
}



### Find all elements:

#' Find all elements related to an ID
#'
#' For a given ID these functions return all IDs of related elements.
#'
#' @details
#'   \code{find_down} finds all elements downwards the hierarchy:
#'
#'   \tabular{rrr}{
#'
#'     node     \tab -> \tab node\cr
#'
#'     way      \tab -> \tab way + node\cr
#'
#'     relation \tab -> \tab relation + way + node\cr
#'
#'   }
#'
#' @param object An \code{\link{osmar}} object
#' @param ids A vector of IDs tagged whether they are \code{node},
#'   \code{way}, or \code{relation}
#'
#' @return A list with the three elements \code{node_ids},
#'   \code{way_ids}, \code{relation_ids}
#'
#' @examples
#'   data("muc", package = "osmar")
#'   o1 <- find(muc, way(tags(k == "highway" & v == "pedestrian")))
#'
#'   find_down(muc, way(o1))
#'   find_up(muc, way(o1))
#'
#' @family finding
#'
#' @export
find_down <- function(object, ids) {
  stopifnot(is_osmar(object))

  handler <- sprintf("find_down_%s", attr(ids, "element"))
  do.call(handler, list(object, as.vector(ids)))
}



find_down_node <- function(object, ids = NULL) {
  stopifnot(is_osmar(object))
  list(node_ids = ids, way_ids = NULL, relation_ids = NULL)
}



find_down_way <- function(object, ids = NULL) {
  ## TODO: check if way id is in object
  stopifnot(is_osmar(object))

  node_ids <- subset_ways(object$ways, ids)$refs$ref
  list(node_ids = node_ids, way_ids = ids, relation_ids = NULL)
}



find_down_relation <- function(object, ids = NULL) {
  ## TODO: check if relation id is in object
  stopifnot(is_osmar(object))

  refs <- subset_relations(object$relations, ids)$refs

  #way_ids <- subset(refs, type == "way")$ref  # CMD check note: no visible binding
  way_ids <- refs[refs$type == "way", ]$ref
  #node_ids <- subset(refs, type == "node")$ref  # CMD check note: no visible binding
  node_ids <- refs[refs$type == "node", ]$ref

  ret <- find_down_way(object, way_ids)
  ret$node_ids <- c(ret$node_ids, node_ids)
  ret$relation_ids <- ids

  ret
}



#' @details
#'   \code{find_up} finds all elements upwards the hierarchy:
#'
#'   \tabular{rrr}{
#'
#'     node     \tab -> \tab node + way + relation\cr
#'
#'     way      \tab -> \tab way + relation\cr
#'
#'     relation \tab -> \tab relation\cr
#'
#'   }
#'
#' @rdname find_down
#'
#' @export
find_up <- function(object, ids) {
  stopifnot(is_osmar(object))

  handler <- sprintf("find_up_%s", attr(ids, "element"))
  do.call(handler, list(object, as.vector(ids)))
}



find_up_node <- function(object, ids = NULL) {
  #way_ids <- subset(object$ways$refs, ref %in% ids)$id  # CMD check note: no visible binding
  way_ids <- object$ways$refs[object$ways$refs$ref %in% ids, ]$id
  #rel_ids <- subset(object$relations$refs, type == "node" & ref %in% ids)$id  # CMD check note: no visible binding
  rel_ids <- object$relations$refs[object$relations$refs$type == "node" &
                                   object$relations$refs$ref %in% ids, ]$id

  list(node_ids = ids, way_ids = way_ids, relation_ids = rel_ids)
}



find_up_way <- function(object, ids = NULL) {
  #rel_ids <- subset(object$relations$refs, type == "way" & ref %in% ids)$id  # CMD check note: no visible binding
  rel_ids <- object$relations$refs[object$relations$refs$type == "way" &
                                   object$relations$refs$ref %in% ids, ]$id
  list(node_ids = NULL, way_ids = ids, relation_ids = rel_ids)
}



find_up_relation <- function(object, ids = NULL) {
  stopifnot(is_osmar(object))
  list(node_ids = NULL, way_ids = NULL, relation_ids = ids)
}



### Find nearest node with given conditions:

#' Find nearest node with given conditions
#'
#' For a given ID, find nearest node (geographical distance) with
#' given conditions.
#'
#' @param object An \code{\link{osmar}} object
#' @param id An node ID
#' @param condition Condition for the element to find; see
#'   \code{\link{find}}
#'
#' @return A node ID or \code{NA}
#'
#' @examples
#'   data("muc", package = "osmar")
#'   id <- find(muc, node(tags(v == "Marienplatz")))[1]
#'
#'   find_nearest_node(muc, id, way(tags(k == "highway" & v == "pedestrian")))
#'
#' @family finding
#'
#' @importFrom geosphere distm
#' @export
find_nearest_node <- function(object, id, condition) {
  stopifnot(is_osmar(object))

  node <- subset_nodes(object$nodes, id)

  element <- attr(condition, "element")

  cand_ids <- find(object, condition)
  cand_ids <- do.call(element, list(cand_ids))
  cand_ids <- find_down(object, cand_ids)
  cand_nodes <- subset_nodes(object$nodes, cand_ids$node_ids)

  dist <- distm(node$attrs[, c("lon", "lat")],
                cand_nodes$attrs[, c("lon", "lat")])

  cand_nodes$attrs[which.min(dist), "id"]
}



#' Binary operators for grep-like functions
#'
#' Binary operators for grep-like functions to use in conditions
#' similar to the "==" operator.
#'
#' @details
#'   \code{%grep%} is currently implemented as \code{grepl(y, x,
#'   ignore.case = TRUE)}.
#'
#' @usage
#'   x %grep% y
#'
#' @rdname binary_grep
#' @export
`%grep%` <- function(x, y) {
  grepl(y, x, ignore.case = TRUE)
}



#' @details
#'   \code{%agrep%} is currently implemented as \code{agrep(y, x,
#'   ignore.case = TRUE)} and converts the index result into a logical
#'   vector.
#'
#' @usage
#'   x %agrep% y
#'
#' @rdname binary_grep
#' @export
`%agrep%` <- function(x, y) {
  w <- agrep(y, x, ignore.case = TRUE)
  r <- rep(FALSE, length(x))
  r[w] <- TRUE
  r
}


