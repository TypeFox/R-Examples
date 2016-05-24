averageFishersZs <- function(zs, ns) {
  if (length(zs) != length(ns)) {
    stop("Vector 'zs' (current length: ", length(zs),
         ") and vector 'ns' (current length: ", length(ns),
         ") must be the same length!");
  }
  return( sum((ns - 3)* zs) / (sum(ns) - 3 * length(ns)));
}