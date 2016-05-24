uLSIF_search_sigma_and_lambda <- function(x, y, centers, sigma_list, lambda_list, verbose) {
  nx <- nrow(x)
  ny <- nrow(y)
  n_min <- min(nx, ny)
  kernel_num <- nrow(centers)

  score_new <- Inf
  sigma_new <- 0
  lambda_new <- 0
  for (sigma in sigma_list) {
    phi_x <- compute_kernel_Gaussian(x, centers, sigma)
    phi_y <- compute_kernel_Gaussian(y, centers, sigma)
    H <- crossprod(phi_y) / ny
    h <- colMeans(phi_x)
    phi_x <- t(phi_x[seq_len(n_min), , drop = FALSE])
    phi_y <- t(phi_y[seq_len(n_min), , drop = FALSE])
    for (lambda in lambda_list) {
      B <- H + diag(lambda * (ny - 1) / ny, nrow = kernel_num, ncol = kernel_num)
      B_inv <- solve(B)
      B_inv_X <- B_inv %*% phi_y
      X_B_inv_X <- phi_y * B_inv_X
      denom <- ones(n_min, value = ny) - ones(kernel_num) %*% X_B_inv_X
      B0 <- B_inv %*% (h %*% ones(n_min)) +
        B_inv_X %*% diag(as.vector((t(h) %*% B_inv_X) / denom))
      B1 <- B_inv %*% phi_x +
        B_inv_X %*% diag(as.vector((ones(kernel_num) %*% (phi_x * B_inv_X)) / denom))
      B2 <- (ny-1) * (nx * B0 - B1) / (ny * (nx - 1))
      B2[B2 < 0] <- 0
      r_y <- t(ones(kernel_num) %*% (phi_y * B2))
      r_x <- t(ones(kernel_num) %*% (phi_x * B2))
      score <- (crossprod(r_y) / 2 - ones(n_min) %*% r_x) / n_min
      if(score < score_new) {
        if(verbose) message(sprintf("  sigma = %.3f, lambda = %.3f, score = %.3f", sigma, lambda, score))
        score_new <- score
        sigma_new <- sigma
        lambda_new <- lambda
      }
    }
  }
  list(sigma = sigma_new, lambda = lambda_new)
}

ones <- function(size, value=1) {
  t(rep(value, size))
}
