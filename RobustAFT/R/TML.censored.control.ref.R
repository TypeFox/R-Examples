TML.censored.control.ref<-function(maxit.sigma=2,tol.sigma=0.0001,maxit.Beta=2,tol.Beta=0.0001,
  Maxit.S=50, tol.S.sigma=0.001, tol.S.Beta=0.001, alg.sigma=1, nitmon = FALSE)
{
  if (!is.numeric(maxit.sigma) || maxit.sigma <= 0)
        stop("value of 'maxit.sigma' must be > 0")
  if (!is.numeric(tol.sigma) || tol.sigma <= 0)
        stop("value of 'tol.sigma' must be > 0")
  if (!is.numeric(maxit.Beta) || maxit.Beta <= 0)
        stop("value of 'maxit.Beta' must be > 0")
  if (!is.numeric(tol.Beta) || tol.Beta <= 0)
        stop("value of 'tol.Beta' must be > 0")
  if (!is.numeric(Maxit.S) || Maxit.S <= 0)
        stop("value of 'Maxit.S' must be > 0")
  if (!is.numeric(tol.S.sigma) || tol.S.sigma <= 0)
        stop("value of 'tol.S.sigma' must be > 0")
  if (!is.numeric(tol.S.Beta) || tol.S.Beta <= 0)
        stop("value of 'tol.S.Beta' must be > 0")
  if (!(alg.sigma %in% c(1, 2)))
        stop("value of 'alg.sigma' must be 1 or 2")
  if (!is.logical(nitmon))
        stop("'nitmon' must be of logical type")

  list(maxit.sigma = maxit.sigma, tol.sigma = tol.sigma, maxit.Beta = maxit.Beta,
    tol.Beta = tol.Beta, Maxit.S = Maxit.S, tol.S.sigma = tol.S.sigma,
    tol.S.Beta = tol.S.Beta, alg.sigma = alg.sigma, nitmon = nitmon)
}