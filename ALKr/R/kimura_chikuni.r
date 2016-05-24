#' @rdname inverse_ALK
#' @export
kimura_chikuni <- function(x, fi1, fi2, threshold = 0.0001, maxiter = 2000,
                          age_classes = colnames(x),
                           length_classes = rownames(x), name = "",
                           description = "") {
  
  invAlk <- calc_invALK(calc_ALK(x), fi1)
  pj2 <- rep(1 / ncol(x), ncol(x))
  criterion <- 1
  iter <- 0
  while (criterion > threshold & iter < maxiter) {
    pj2.old <- pj2
    alk <- t(t(invAlk) * pj2) / rowSums(t(t(invAlk) * pj2))
    nij <- fi2 * alk
    pj2 <- colSums(nij) / sum(nij)
    criterion <- sum(abs(pj2 - pj2.old))
    iter <- iter + 1
  }
  
  new("ALKr", alk = calc_ALK(nij),
      N = nij,
      method = "Kimura & Chikuni",
      parameters = list(
        ConvergenceThreshold = threshold,
        Iterations = iter,
        Converged = iter < maxiter),
      name = name,
      description = description
  )
}

