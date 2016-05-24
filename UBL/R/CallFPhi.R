
# ======================================================================
# This function calls the Fortran code for defining and evaluating the 
# relevance function "phi".
# P.Branco Apr 2016
# ======================================================================

phi <- function(y, control.parms) {
  n <- length(y)
  Charmeth <- control.parms[[1]]
  meth <- ifelse(Charmeth == "extremes", 0, 1) 
  npts <- control.parms[[2]]
  lparms <- length(control.parms[[3]])
  phiParms <- control.parms[[3]]
  yPhi <- rep(0.0, times=n)
  ydPhi <- rep(0.0, times=n)
  yddPhi <- rep(0.0, times=n)  

  storage.mode(n) <- "integer"
  storage.mode(y) <- "double"
  storage.mode(meth) <- "integer"
  storage.mode(npts) <- "integer"
  storage.mode(lparms) <- "integer"
  storage.mode(phiParms) <- "double"
  storage.mode(yPhi) <- "double"
  storage.mode(ydPhi) <- "double"
  storage.mode(yddPhi) <- "double"
  

  res <- .Fortran("rtophi",
            n = n, # nr of points
            y = y,  # tgt values
            method = meth, # coded method (0:extremes; 1:range)
            npts = npts, #the nr of points in the relevance matrix provided
            lparms = lparms,  # length of phiParms
            phiParms = phiParms, # matrix info
            yPhi = yPhi,# output
            ydPhi = ydPhi,# output not used
            yddPhi = yddPhi# output not used
            )
  
res <- res$yPhi

}
