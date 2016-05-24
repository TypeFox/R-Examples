##' This is an implmentation of the FIRE algorithm for structural
##' relaxation put forward by Bitzek et al. (2006)
##' @title The FIRE algorithm
##' @param r Initial locations of particles
##' @param force Force function
##' @param restraint Restraint function
##' @param m Masses of points
##' @param dt Initial time step
##' @param maxmove Maximum distance to move in any time step
##' @param dtmax Maxiumum time step
##' @param Nmin Number of steps after which to start increasing \code{dt}
##' @param finc Fractional increase in \code{dt} per time step
##' @param fdec Fractional decrease in \code{dt} after a stop
##' @param astart Starting value of \code{a} after a stop
##' @param fa Fraction of \code{a} to retain after each step
##' @param a Initial value of \code{a}
##' @param nstep Maxiumum number of steps
##' @param tol Tolerance - if RMS force is below this value, stop and
##' report convergence
##' @param verbose If \code{TRUE} report progress verbosely
##' @return List containing \code{x}, the positions of the points,
##' \code{conv}, which is 0 if convergence as occured and 1 otherwise,
##' and \code{frms}, the root mean square of the forces on the
##' particles.
##' @references Bitzek, E., Koskinen, P., G\"{a}hler, F., Moseler, M.,
##' and Gumbsch, P. (2006). Structural relaxation made
##' simple. Phys. Rev. Lett., 97:170201.
##' @author David Sterratt
fire <- function(r, force, restraint, m=1, dt=0.1, maxmove=1E2, dtmax=1,
                 Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99, a=0.1,
                 nstep=100, tol=0.00001, verbose=FALSE) {
  Nsteps <- 0
  conv <- 1
  # Initialise velocity
  v <- 0*r                              
  dt <- min(dt, dtmax)

  ## Counters for number of stops and hits of dtmax and maxmove
  nstop <- 0
  ndtmax <- 0
  nmaxmove <- 0

  for (i in 1:nstep) {
    f <- force(r)
    frad2 <- dot(f, r/vecnorm(r), 2)^2
    ## print(frad2)
    ftan2 <- dot(f, f, 2) - frad2
    ##print(ftan2)
    frms <- sqrt(mean(dot(f, f, 2)))
    ftanrms <- sqrt(mean(ftan2))
    if (ftanrms < tol) {
      conv <- 0
      break;
    }
    vf <- dot(f, v, 2)
    ## v[vf>0,] <- 0
    if (sum(vf) > 0) {
      v <- (1 - a)*v + a*f/vecnorm(f)*vecnorm(v)
      if (Nsteps > Nmin) {
        dt <- min(dt*finc, dtmax)
        if (dt==dtmax) ndtmax <- ndtmax + 1 
        a <- a*fa
      }
      Nsteps <- Nsteps + 1
    } else {
      nstop <- nstop +1
      v <- 0
      a <- astart
      dt <- dt*fdec
      Nsteps <- 0
    }
    v <- v + dt*f/m
    dr <- dt*v
    normdr <- sqrt(sum(dr*dr))
    if (normdr > maxmove) {
      if (verbose) nmaxmove <- nmaxmove + 1
      dr <- maxmove*dr/normdr
    }
    ## if (!is.null(mm)) {
    ##   normdr <- vecnorm(dr)
    ##   dr[normdr>mm,] <- dr[normdr>mm,]/normdr[normdr>mm]*mm[normdr>mm]
    ## }
    r <- r + dr
    r <- restraint(r)
  }
  if (verbose) {
    message(paste("FIRE: ", nstop, " stops. ",
                  nmaxmove, "hits of maxmove. ",
                  ndtmax, "hits of dtmax."))
    message("Frms = ", frms, "; Ftanrms = ", ftanrms)
  }
  return(list(x=r, conv=conv, frms=frms))
}

