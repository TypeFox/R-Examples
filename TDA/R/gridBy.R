gridBy <-
function(lim = c(0, 1), by = (lim[2] - lim[1]) / 10) {
  gridlist <- list()
  for (idx in 1:(length(lim)/2)) {
    if (length(by) == 1) {
      gridlist[[idx]] <- seq(lim[2 * idx - 1], lim[2 * idx], by = by)
    } else {
      gridlist[[idx]] <- seq(lim[2 * idx - 1], lim[2 * idx], by = by[idx])
    }
  }
  grid <- expand.grid(gridlist)
  colnames(grid) <- NULL
  dim <- sapply(gridlist, length)
  out <- list("dim" = dim, "lim" = lim, "grid" = grid)
  return(out)
}