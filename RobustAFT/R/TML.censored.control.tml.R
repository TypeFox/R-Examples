TML.censored.control.tml<-function(maxit.sigma=20,tol.sigma=0.0001,maxit.Beta=20,tol.Beta=0.0001,
           Maxit.TML=50, tol.TML.sigma=0.001, tol.TML.Beta=0.001, alg.sigma=1, nitmon = FALSE)
{
  if (!is.numeric(maxit.sigma) || maxit.sigma <= 0)
        stop("value of 'maxit.sigma' must be > 0")
  if (!is.numeric(tol.sigma) || tol.sigma <= 0)
        stop("value of 'tol.sigma' must be > 0")
  if (!is.numeric(maxit.Beta) || maxit.Beta <= 0)
        stop("value of 'maxit.Beta' must be > 0")
  if (!is.numeric(tol.Beta) || tol.Beta <= 0)
        stop("value of 'tol.Beta' must be > 0")
  if (!is.numeric(Maxit.TML) || Maxit.TML <= 0)
        stop("value of 'Maxit.TML' must be > 0")
  if (!is.numeric(tol.TML.sigma) || tol.TML.sigma <= 0)
        stop("value of 'tol.TML.sigma' must be > 0")
  if (!is.numeric(tol.TML.Beta) || tol.TML.Beta <= 0)
        stop("value of 'tol.S.Beta' must be > 0")
  if (!(alg.sigma %in% c(1, 2)))
        stop("value of 'alg.sigma' must be > 0")
  if (!is.logical(nitmon))
        stop("'nitmon' must be of logical type")

  list(maxit.sigma = maxit.sigma, tol.sigma = tol.sigma, maxit.Beta = maxit.Beta,
    tol.Beta = tol.Beta, Maxit.TML = Maxit.TML, tol.TML.sigma = tol.TML.sigma,
    tol.TML.Beta = tol.TML.Beta, alg.sigma = alg.sigma, nitmon = nitmon)
}