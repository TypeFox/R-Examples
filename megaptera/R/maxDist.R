maxDist <- function(s, tolerance = 0) {
  d <- dist.dna(s, model = "N", pairwise.deletion = TRUE,
                as.matrix = TRUE)
  diag(d) <- NA
  min.dist <- min(d, na.rm = TRUE)
  if ( min.dist > tolerance ) out <- min.dist
  else {
    diag(d) <- 0
    id <- which.min(apply(d, 2, sum, na.rm = TRUE))
    id <- d[, id] <= tolerance
    names(id) <- NULL
    out <- if ( any(id) ) s[id, ]
  }
  out
}

getMaxDist <- function(s, tolerance = 0) {
  d <- dist.dna(s, model = "N", pairwise.deletion = TRUE,
                as.matrix = TRUE)
  max(d, na.rm = TRUE)
}