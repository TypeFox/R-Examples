#' @rdname inverse_ALK
#' @export
hoenig_heisey <- function(x, fi1, fi2, threshold = 0.0001, maxiter = 2000,
                         age_classes = colnames(x),
                          length_classes = rownames(x), name = "",
                          description = "") {
  
  nij1 <- fi1 * calc_ALK(x)
  nij2 <- fi2 * calc_ALK(x)
  pj2 <- colSums(nij2) / sum(nij2)
  
  criterion <- 1
  iter <- 0
  
  while(criterion > threshold & iter < maxiter) {
    pj2.old <- pj2
    ialk <- sweep(nij1 + nij2, 2, colSums(nij1 + nij2), "/")
    ialk[is.nan(ialk)] <- 0
    numer <- sweep(ialk, 2, colSums(nij2), "*")
    denom <- rowSums(numer)
    denom[denom == 0] <- 1
    nij2 <- fi2 * sweep(numer, 1, denom, "/")
    pj2 <- colSums(nij2) / sum(nij2)
    criterion <- sum(abs(pj2 - pj2.old))
    iter <- iter + 1
  }
  
  new("ALKr", alk = calc_ALK(nij2),
      N = nij2,
      method = "Hoenig & Heisey",
      parameters = list(
        ConvergenceThreshold = threshold,
        Iterations = iter,
        Converged = iter < maxiter),
      name = name,
      description = description
  )
}
