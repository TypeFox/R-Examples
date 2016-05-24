#' @include as-osmar.R
#' @include osmar-plotting.R
{}



#' Convert osmar object to igraph
#'
#' Convert an osmar object to an igraph (see
#' igraph-package).
#'
#' @param obj An \code{\link{osmar}} object
#'
#' @return An igraph-package \code{graph} object
#' 
#' @examples
#' file <- system.file("extdata", "kaufstr.xml", package = "osmar")
#' raw <- readLines(file)
#' kaufstr <- as_osmar(xmlParse(raw))
#' kaufstrGraph <- as_igraph(kaufstr)
#'
#' @export
as_igraph <- function(obj) {
  stopifnot(is_osmar(obj))
  stopifnot(require("igraph"))

  dat <- merge_ways_nodes(obj$ways[[3]], obj$nodes[[1]])
  dat <- split(dat, dat$id)
  dat <- dat[sapply(dat, nrow) >= 2]

  edges <- lapply(dat,
                  function(x) {
                    n <- nrow(x)
                    from <- 1:(n-1)
                    to <- 2:n

                    weights <- distHaversine(x[from, c("lon", "lat")],
                                             x[to, c("lon", "lat")])

                    cbind(from_node_id = x[from, "ref"],
                          to_node_id = x[to, "ref"],
                          way_id = x[1, "id"],
                          weights = weights)
                  })
  edges <- do.call(rbind, edges)

  weights <- edges[, "weights"]
  names <- edges[, "way_id"]
  edges <- cbind(as.character(edges[, "from_node_id"]),
                 as.character(edges[, "to_node_id"]))

  graph <- graph.edgelist(edges)
  E(graph)$weight <- weights
  E(graph)$name <- names

  graph
}


## as_igraph <- function(obj) {
##   way_nodes <- split(obj$ways[[3]]$ref, obj$ways[[3]]$id)
##   way_nodes <- ways[sapply(way_nodes, length) >= 2]

##   edges2 <- lapply(way_nodes,
##          function(x) {
##            dat <- subset_nodes(obj$nodes, x)$attrs

##            n <- length(x)
##            from <- 1:(n-1)
##            to <- 2:n

##            weights <- distHaversine(dat[from, c("lon", "lat")],
##                                     dat[to, c("lon", "lat")])

##            cbind(from_node_id = dat[from, "id"],
##                  to_node_id = dat[to, "id"],
##                  weights = weights)
##          })
##   edges2 <- do.call(rbind, edges2)
##   edges2 <- na.omit(edges2)

##   weights <- edges[, "weights"]
##   edges <- cbind(as.character(edges[, "from_node_id"]),
##                  as.character(edges[, "to_node_id"]))

##   graph <- graph.edgelist(edges)
##   E(graph)$weight <- weights

##   graph

## }
