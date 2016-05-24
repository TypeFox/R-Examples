# helper function counting peaks, takes threshold
get_k_n <- function(x, threshold = 1) {
  #   if (is.null(threshold)) {
  #     x <- x[x != 0]
  #     threshold <- median(x)
  #   }
  peaks <- findpeaks(x, threshold = threshold)
  pos <- nrow(peaks)
  if (is.null(pos))
    pos <- 0
  # list(k = pos, n = nrow(findpeaks(x, threshold = min(x))), thr = threshold)
  pos
}