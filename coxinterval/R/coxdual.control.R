coxdual.control <- function(eps = 1e-7, iter.max = 50000, coef.typ = 1,
                            coef.max = 10, sieve = TRUE, sieve.const = 1,
                            sieve.rate = 1/3, risk.min = 1, data = FALSE)
{
  if (eps <= .Machine$double.eps)
    stop("Invalid epsilon. Choose a small value > ", .Machine$double.eps, ".")
  if (iter.max < 0)
    stop("Invalid maximum iterations. Choose a large positive integer.")
  if (coef.typ < eps)
    stop("Invalid coefficient magnitude. Choose a positive value.")
  if (coef.max <= coef.typ)
    stop("Invalid maximum coefficient size. Choose a value > ", coef.typ, ".")
  if (any(sieve.const < eps))
    stop("Invalid sieve constant. Choose a positive value.")
  if (sieve & length(sieve.const != 3))
    sieve.const <- rep(sieve.const[1], 3)
  if (sieve & (sieve.rate <= 1/8 | sieve.rate >= 1/2))
    stop("Invalid sieve rate. Choose a value in (1/8, 1/2).")
  if (risk.min < 1 | risk.min > 2)
    stop("Minimum risk set should be 1 or 2.")
  list(eps = eps, iter.max = iter.max, coef.typ = coef.typ, coef.max = coef.max,
       sieve = sieve, sieve.const = sieve.const, sieve.rate = sieve.rate,
       risk.min = risk.min, data = data)
}
