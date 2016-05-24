KLIEP_search_sigma_list <- function(x, y, centers, sigma_list, fold, verbose = TRUE) {
  nx <- nrow(x)

  score <- -Inf

  for(sigma_new in sigma_list) {
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
    if (score_new > score) {
      score <- score_new
      sigma <- sigma_new
      if(verbose) message(sprintf("  sigma = %.5f, score = %f", sigma, score))
    }
  }
  sigma
}
