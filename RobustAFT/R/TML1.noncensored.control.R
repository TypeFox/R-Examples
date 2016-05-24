TML1.noncensored.control <- function(iv = 1, gam = 0.1, maxit = 200, tol = 1e-04)
{
  if (!is.numeric(iv) || !(iv %in% c(0, 1)))
        stop("value of 'iv' must be 0 or 1")
  if (!is.numeric(gam) || gam <= 0 || gam > 1)
        stop("value of 'gam' must be between 0 and 1")
  if (!is.numeric(maxit) || maxit < 2)
        stop("value of 'maxit' must be an integer > 1")
  if (!is.numeric(tol) || tol <= 0)
        stop("'tol' must be > 0 ")
  list(iv = iv, gam = gam, maxit=maxit, tol=tol)
}


