TML.noncensored.control<-function(iv = 1, nrep = 0, gam = 0.1, nitmon = FALSE, maxit = 200, tol = 1e-04, fastS = FALSE, seed=1313)
{
  if (!is.numeric(iv) || !(iv %in% c(0, 1)))
        stop("value of 'iv' must be 0 or 1")
  if (!is.numeric(nrep) || nrep < 0)
        stop("value of 'nrep' must be 0 or higher")
  if (!is.numeric(gam) || gam <= 0 || gam > 1)
        stop("value of 'gam' must be between 0 and 1")
  if (!is.logical(nitmon))
        stop("'nitmon' must be of logical type")
  if (!is.numeric(maxit) || maxit < 2)
        stop("value of 'maxit' must be an integer > 1")
  if (!is.numeric(tol) || tol <= 0)
        stop("'tol' must be > 0 ")
 if (!is.logical(fastS))
        stop("'fastS' must be of logical type")
  if (!is.numeric(seed) || seed < 1)
        stop("'seed' must be an integer > 0r")  
  list(iv = iv, nrep=nrep, gam = gam, nitmon=nitmon, maxit=maxit, tol=tol, fastS=fastS, seed=seed)
}


