simulate.transectHolder <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  distances <- c()
  while(length(distances) < nsim) {
    unfiltered <- do.call(object$rng,
                          as.list(c(10*nsim, object$parameters)))
    filter <- runif(10*nsim, 0, max(unfiltered))
    distances <- c(distances, unfiltered[unfiltered > filter])
    }
  distances <- distances[1:nsim]
  angles <- runif(nsim, 0, 2*pi)
  return(data.frame(distances = distances,
                    angles = angles,
                    x = cos(angles) * distances,
                    y = sin(angles) * distances))
}
