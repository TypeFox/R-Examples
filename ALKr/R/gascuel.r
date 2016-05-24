#' @rdname inverse_ALK
#' @useDynLib ALKr
#' @export
gascuel <- function(x, fi1, fi2, initial_values,
                    threshold = 0.0001, maxiter = 2000,
                    age_classes = colnames(x), length_classes = rownames(x),
                    name = "", description = "") {
  
  nij <- fi1 * calc_ALK(x)
  li <- as.numeric(length_classes)
  lj <- colSums(nij * li) / colSums(nij)
  
  pi_ <- fi1 / sum(fi1)
  
  optimal <- optim(initial_values, fn = optimGascuel, lj = lj, li = li,
                   pi_ = pi_, threshold = threshold, maxiter = maxiter)
  
  if (optimal$convergence == 1)
    warning("Parameter optimization exceeded maxiter")
  if (optimal$convergence == 10)
    warning("Degeneracy of the Nelder-Mead simplex optimization")
  
  result <- calc_ALK(finalGascuel(optimal$par, lj = lj, li = li, pi_ = pi_,
                                  threshold = threshold, maxiter = maxiter) * sum(fi2))
  
  rownames(result) <- rownames(nij)
  colnames(result) <- colnames(nij)
  
  new("ALKr",
      alk = result,
      N = nij,
      method = "Gascuel",
      parameters = list(
        ConvergenceThreshold = threshold,
        alpha = optimal$par[1],
        beta = optimal$par[2],
        gamma = optimal$par[3],
        Converged = optimal$convergence == 0),
      name = name,
      description = description
  )
}
