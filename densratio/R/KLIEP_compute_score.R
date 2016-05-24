KLIEP_compute_score <- function(phi, alpha) {
  mean(log(phi %*% alpha))
}
