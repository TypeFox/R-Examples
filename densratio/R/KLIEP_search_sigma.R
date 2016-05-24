KLIEP_search_sigma <- function(x, y, centers, fold, verbose = TRUE) {
  nx <- nrow(x)

  sigma <- 10
  score <- -Inf

  for (digit_position in 0:-5) {
    for (i in 1:9) {
      sigma_new <- sigma - 10 ^ digit_position
      cv_split <- sample(nx) %% fold + 1

      phi_x <- compute_kernel_Gaussian(x, centers, sigma_new)
      phi_y <- compute_kernel_Gaussian(y, centers, sigma_new)
      mean_phi_y <- matrix(colMeans(phi_y))

      scores <- numeric(fold)
      for (i in seq_len(fold)) {
        alpha <- KLIEP_optimize_alpha(phi_x[cv_split != i, ], mean_phi_y = mean_phi_y)
        scores[i] <- KLIEP_compute_score(phi_x[cv_split == i, ], alpha)
      }
      score_new <- mean(scores)
      if (score_new <= score) {
        break
      }
      score <- score_new
      sigma <- sigma_new
      if(verbose) message(sprintf("  sigma = %.5f, score = %f", sigma, score))
    }
  }
  sigma
}
