ordered.clusters <- function(id) {
  d <- which(duplicated(id)) # NB: d > 1, poss. empty
  all(id[d]==id[d-1])
}
