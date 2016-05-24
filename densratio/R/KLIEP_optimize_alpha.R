KLIEP_optimize_alpha <- function(phi_x, phi_y, mean_phi_y) {
  A <- phi_x
  if(missing(phi_y)) {
    if(is.matrix(mean_phi_y)) {
      b <- mean_phi_y
    } else {
      b <- matrix(mean_phi_y)
    }
  } else {
    b <- matrix(colMeans(phi_y))
  }
  c <- b / crossprod(b)[1,1]

  kernel_num <- ncol(phi_x)
  max_iteration <- 100

  alpha <- matrix(rep(1, kernel_num))
  alpha <- compute_next_alpha(alpha, b, c)
  score <- KLIEP_compute_score(phi_x, alpha)

  eplilon_list <- 10 ^ (3:-3)
  for (epsilon in eplilon_list) {
    for (i in seq_len(max_iteration)) {
      alpha_new <- alpha + (epsilon * t(A)) %*% (1 / (A %*% alpha))
      alpha_new <- compute_next_alpha(alpha_new, b, c)
      score_new <- KLIEP_compute_score(phi_x, alpha_new)
      if(score_new <= score) {
        break
      }
      alpha <- alpha_new
      score <- score_new
    }
  }
  alpha
}

compute_next_alpha <- function(alpha, b, c = b / crossprod(b)[1,1]) {
  alpha <- alpha + (1 - (t(b) %*% alpha)[1,1]) * c
  alpha[alpha < 0, ] <- 0
  alpha <- alpha / (t(b) %*% alpha)[1,1]
  alpha
}

