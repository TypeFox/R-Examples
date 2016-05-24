heckitrob.control <-
function(acc = 1e-04, test.acc = "coef", maxit = 50, maxitO = 50, weights.x1 = c("none", "hat", "robCov", "covMcd"),
                            weights.x2 = c("none", "hat", "robCov", "covMcd"), tcc = 1.345, t.c = 1.345)
{
  if (!is.numeric(acc) || acc <= 0) 
    stop("value of acc must be > 0")
  if (test.acc != "coef") 
    stop("Only 'test.acc = \"coef\"' is currently implemented")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(maxitO) || maxitO <= 0) 
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(tcc) || tcc <= 0) 
    stop("value of the tuning constant c (tcc) must be > 0")
  if (!is.numeric(t.c) || t.c <= 0) 
    stop("value of the tuning constant c (t.c) must be > 0")
  if (!is.character(weights.x1)) 
    stop("choose the implemented method of the weight function")
  if (!is.character(weights.x2)) 
    stop("choose the implemented method of the weight function")
  list(acc = acc, test.acc = test.acc, maxit = maxit, maxitO = maxitO, weights.x1 = weights.x1, weights.x2 = weights.x2[1], tcc = tcc, t.c = t.c)
}
