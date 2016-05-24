TML.censored.control.S <-function(N = 100, q=6, sigma0 = 1, MAXIT = 100, TOL = 0.01, seed = 153)
{
  if (!is.numeric(N) || N < 21)
        stop("value of 'N' must be > 20")
  if (!is.numeric(q) || q < 1)
        stop("value of 'q' must be > 0")
  if (!is.numeric(sigma0) || sigma0 <= 0 )
        stop("value of 'sigma0' must be > 0")
  if (!is.numeric(MAXIT) || MAXIT < 2)
        stop("value of 'MAXIT' must be an integer > 1")
  if (!is.numeric(TOL) || TOL <= 0)
        stop("'TOL' must be > 0 ")
  if (!is.numeric(seed) || seed < 1)
        stop("'seed' must be an integer > 0")  
  list(N = N, q=q, sigma0 = sigma0, MAXIT=MAXIT, TOL=TOL, seed=seed)
}

