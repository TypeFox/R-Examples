jackknife.smacofB <- function(object, eps = 1e-6, itmax = 100, verbose = FALSE) 
{
  ## object... object of class smacofB (from smacofSym, smacofConstraint)
    if (class(object)[1] != "smacofB") stop("Jackknife is currenlty implemented for objects of class smacofB from smacofSym() only! \n")
    if (object$model == "SMACOF constraint") stop("Jackknife is currenlty implemented for smacofSym() objects only! \n")
    
    delta <- as.matrix(object$delta)
    n <- nrow(delta)
    x0 <- object$conf
    ndim <- object$ndim
   type <- object$type
    
    #x0 <- smacofSym (delta, ndim = ndim, metric = metric) $ conf
    xx <- smacofDeleteOne(delta, ndim = ndim, type = type)
    kk <- array(rep(diag (ndim), n), c(ndim, ndim, n))
    cc <- matrix(0, n, ndim)
    bb <- matrix(0, n, ndim)
    yy <- xx
    oloss <- Inf
    itel <- 1
    repeat {
        y0 <- matrix (0, n, ndim)
        for (i in 1 : n) {
            y0 <- y0 + xx[, , i] %*% kk[, , i]
            }
        y0 <- ((n - 1) * y0) / (n * (n - 2))
        for (i in 1: n) {
            zz <- matrix (0, n, ndim)
            for (j in i : n) {
                zz <- zz + xx[, , j] %*% kk[, , j]
                }
            xz <- crossprod ( xx[, , i], zz)
            kk[, , i] <- procrustus (xz)
            }
        nloss <- 0
        for (i in 1 : n) {
  
 
            yy[, , i] <- xx[, , i] %*% kk [, , i] 
            yy[i, , i] <- n * y0[i, ] / (n - 1)
            yy[, , i] <- yy[, , i] - outer (rep (1, n), y0[i, ] / (n - 1))
            nloss <- nloss + sum ( (y0 - yy[, , i]) ^ 2)
            }
        if (verbose) {
  	     cat("Iteration: ", formatC (itel, digits=3, width=3),
			 "Old Loss: ", formatC (oloss, digits=10, width=15, format="f"),
			 "New Loss: ", formatC (nloss, digits=10, width=15, format="f"),
			 "\n")
			 }
	    if (((oloss - nloss) < eps) || (itel == itmax)) {
	        break()
	        }
	    itel <- itel + 1
	    oloss <- nloss
	  }
	 x0 <- x0 %*% procrustus(crossprod (x0, y0))
   
   if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
	
  ## stability measure
  stab.num <- sum(apply(yy, 3, function(yyi) (norm(yyi-y0))^2))  ## numerator stab
  stab.denom <- sum(apply(yy, 3, function(yyi) (norm(yyi))^2))   ## denominator stab
  stab <- 1 - stab.num/stab.denom

  ## cross vallidity
  cross.num <- n*(norm(x0-y0))^2
  cross <- 1 - cross.num/stab.denom
  
  ## dispersion around x0
  disp <- 2 - (stab + cross)
  
	result <- list(smacof.conf = x0, jackknife.conf = yy, comparison.conf = y0, stab = stab, cross = cross, disp = disp, niter = itel, loss = nloss, nobj = n, call = match.call())
  class(result) <- "smacofJK"
  result
}
    