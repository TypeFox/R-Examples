#' Example clusterfly object created with olives data
#'
#' @keywords dataset
#' @export
olive_example <- function() {
  ol <- clusterfly(olives[, -(1:3)], olives[, 2:3])
  ol <- cfly_cluster(ol, kmeans, 3, name="kmeans")
  ol <- cfly_cluster(ol, kmeans, 4, name="k4-1")
  ol <- cfly_cluster(ol, kmeans, 4, name="k4-2")
  ol <- cfly_cluster(ol, kmeans, 4, name="k4-3")
  ol[["Region"]] <- olives$Region
  ol[["Area"]] <- olives$Area

  ol
}
