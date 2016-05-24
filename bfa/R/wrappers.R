# @nord
.updateScores <- function (Z_, A_, F_) 
  .Call("updateScoresC", Z_, A_, F_, PACKAGE = "bfa")

# @nord
.updateRho <- function (rho_, A_, rhoa_, rhob_) 
  .Call("updateRho", rho_, A_, rhoa_, rhob_, PACKAGE = "bfa")

# @nord
#.updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_, SEXP A_restrict_, SEXP pnz_ )
.updateSparseLoadings <- function (Z_, A_, F_, tauinv_, rho_, A_restrict_, pnz_) 
  .Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, A_restrict_, pnz_, PACKAGE = "bfa")

#SEXP updateZ( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ ){

# @nord
.updateZ <- function (Z_, Ra_, maxes_, argsorts_, A_, F_) {
  .Call("updateZ", Z_, Ra_, maxes_, argsorts_, A_, F_,  PACKAGE = "bfa")
  return()
}

# @nord
.updateZcut <- function (Z_, Ra_, maxes_, argsorts_, A_, F_) {
  .Call("updateZcut", Z_, Ra_, maxes_, argsorts_, A_, F_,  PACKAGE = "bfa")
  return()
}

# @nord
#SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_)
.MCMCstep <- function (Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, A_restrict_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings, more_args)
  .Call("MCMCstep", Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, A_restrict_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings, more_args, PACKAGE = "bfa")
