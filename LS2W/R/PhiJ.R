`PhiJ` <-
function(J, filter.number = 10., family = "DaubLeAsymm", tol = 1e-100, OPLENGTH
 = 2000., verbose = FALSE)
{
#
# Program to compute the discrete autocorrelation wavelets (1-D)
# This is essentially the same as in WaveThresh, except that I have adopted a 
# different naming convention to distinguish between 1D and 2D autocorrelation
# wavelets (see Psi1Dname or Psi2Dname)
#
if(verbose) {
   cat("Computing PhiJ\n")
   now <- proc.time()[1.:2.]
}
if(J >= 0.)
   stop("J must be a negative integer!")
if(J - round(J) != 0.)
   stop("J must be an integer")
#
# See if matrix already exists
#
Phiorig <- #
Phi1Dname(J = J, filter.number = filter.number, family = family)
{
#
# If it doesn't, go create the father autocorrelation wavelets.
# This uses the routine PhiJ contained within orig.c
#
   if(exists(Phiorig)) {
      if(verbose)
         cat("Returning precomputed version \n")
      if(verbose) {
         speed <- proc.time()[1.:2.] - now
      cat("Took ", sum(speed), " seconds \n")
      }
      return(get(Phiorig,envir=DWEnv)) 
   }
   H <- filter.select(filter.number = filter.number, family = 
   family)$H
   wout <- rep(0., OPLENGTH)
   rlvec <- rep(0.,  - J)
   error <- 0.
   answer <- .C("PhiJ",
   J = as.integer( - J),
   H = as.double(H),
   LengthH = as.integer(length(H)),
   tol = as.double(tol),
   wout = as.double(wout),
   lwout = as.integer(length(wout)),  
   rlvec = as.integer(rlvec),
   error = as.integer(error),
   PACKAGE = "LS2W")
   if(answer$error != 0.) {
      if(answer$error == 160.)
         cat("Increase OPLENGTH to be larger than ", answer$lwout,"\n")
      stop(paste("Error code was ", answer$error))
   }
   if(verbose) {
      speed <- proc.time()[1.:2.] - now
      cat("Took ", sum(speed), " seconds \n")
   }
   m <- vector("list",  - J)
   lj <- c(0., cumsum(2. * answer$rlvec - 1.))
   for(j in 1.:( - J))
      m[[j]] <- answer$wout[(lj[j] + 1.):lj[j + 1.]]
   assign(Phiorig, m, envir=DWEnv)
   m
   }
}

