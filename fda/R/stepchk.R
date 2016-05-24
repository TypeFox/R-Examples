stepchk <- function(oldstep, cvec, deltac, limwrd, ind,
                    climit=50*c(-rep(1,ncvec), rep(1,ncvec)),
                    active=1:ncvec, dbgwrd){
# last changed 2007 or 2008 by Spencer Graves
#  check the step size to keep parameters within boundaries
        ncvec   <- length(deltac)
        bot     <- climit[1:ncvec]
        top     <- climit[ncvec+(1:ncvec)]
        limflg  <- FALSE
        newstep <- oldstep
        stepi   <- oldstep*deltac
        stepmin <- min(stepi)
        index   <- stepi[active] == stepmin
#  ensure that step does not go beyond lower limit on parameters
        if (any(stepi[index] < bot[index]-cvec[index]) &
            any(deltac[index] != 0) )  {
          anew <- min((bot[index]-cvec[index])/deltac[index])
          if (dbgwrd) {
            print("Lower limit reached ... new step:")
            cat(c(stepi, round(c(oldstep, anew),4)),"\n")
            cat(round(cvec + anew*deltac,4),"\n")
          }
          newstep <- anew
          limflg <- TRUE
        }
#  ensure that step does not go beyond upper limit on parameters
        stepi   <- oldstep*deltac
        stepmax <- max(stepi)
        index   <- stepi[active] == stepmax
        if (any(stepi[index] > top[index]-cvec[index]) &
            any(deltac[index] != 0) ) {
          anew <- min((top[index]-cvec[index])/deltac[index])
          if (dbgwrd) {
            print("Upper limit reached ... new step:")
            cat(c(stepi, round(c(oldstep, anew),4)),"\n")
          }
          newstep <- anew
          limflg <- TRUE
        }
#  check whether lower limit has been reached twice in a row
        if (limflg) {
          if (limwrd) ind <- 1 else limwrd <- TRUE
        } else limwrd <- FALSE
  return(list(newstep, ind, limwrd))
}
