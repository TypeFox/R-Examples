`trapTransect` <-
function(distances, seed.counts, trap.area = 0.0001) {
  if (length(distances) == length(seed.counts)) {
    n <- length(distances)
  } else {
    stop("Lengths of distances and counts differ.")
  }
  if ((length(trap.area) > 1)) stop("Ambiguous trap area.")
  trapTransect <- list(distances = distances,
                    seed.counts = seed.counts,
                    trap.area = trap.area)
  class(trapTransect) <- c("trapTransect","list")
  return(trapTransect)
}

