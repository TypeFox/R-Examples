.graythresh <- function(I) {
  I <- as.vector(I)
  if (min(I) < 0 | max(I) >1) {
    stop("Data needs to be between 0 and 1")
  }
  I <- I * 256
  
  num_bins <- 256
  counts <- hist(I,num_bins, plot = FALSE)$counts
  
  # Variables names are chosen to be similar to the formulas in the Otsu paper.
  p <-  counts / sum(counts)
  omega <- cumsum(p)
  mu <- cumsum(p * (1:num_bins))
  mu_t <- mu[length(mu)]

  sigma_b_squared <- (mu_t * omega - mu)^2 / (omega * (1 - omega))
  sigma_b_squared[is.na(sigma_b_squared)] <- 0
  maxval <- max(sigma_b_squared)
  isfinite_maxval <- is.finite(maxval)
  if (isfinite_maxval) {
      idx <- mean(which(sigma_b_squared == maxval))
      #Normalize the threshold to the range [0, 1].
      level <- (idx - 1) / (num_bins - 1)
  } else {
      level <- 0.0
  }
}
