phfit.3mom <- function(m1, m2, m3, method = c("Osogami06", "Bobbio05"),
	max.phase = 50, epsilon = sqrt(.Machine$double.eps)) {
  method <- match.arg(method)
  if (method == "Osogami06") {
    res <- matching3PH(c(m1, m2, m3), epsilon=epsilon, max.phase=max.phase)
    ph(alpha=res$tau, Q=res$T, xi=res$xi)
  } else {
    mm.bobbio05(m1, m2, m3)
  }
}
