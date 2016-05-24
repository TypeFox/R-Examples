### set parameters controlling the model fit
coxaalen.control <- function(eps = 1e-7, eps.norm = c("max", "grad"),
                               iter.max = 5000, armijo = 1/3, var.coef = TRUE,
                               coef.typ = 1, coef.max = 10, trace = FALSE,
                               thread.max = 1, data = FALSE)
{
  eps.norm <- match.arg(eps.norm)
  if (eps <= .Machine$double.eps)
    stop("Invalid epsilon. Choose a small value > ", .Machine$double.eps, ".")
  if (!is.element(eps.norm, c("max", "grad")))
    stop(paste("Unknown stopping rule norm", eps.norm))
  if (iter.max < 1)
    stop("Invalid maximum iterations. Choose a large positive integer.")
  if (var.coef & coef.typ < eps)
    stop("Invalid coefficient magnitude. Choose a positive value.")
  if (var.coef & coef.max <= coef.typ)
    stop("Invalid maximum coefficient size. Choose a value > ", coef.typ, ".")
  if (armijo < eps | armijo >= 1/2)
    stop("Invalid scale for Armijo's rule. Choose a value in (0, 1/2).")
  list(eps = eps, eps.norm = eps.norm, iter.max = iter.max, armijo = armijo,
       var.coef = var.coef, coef.typ = coef.typ, coef.max = coef.max,
       trace = trace, thread.max = thread.max, data = data)
}
