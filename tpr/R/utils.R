tpr.control <- function(conlist) {
  control <- list(maxit = 25, tol = 0.0001, smooth = 0, intsmooth = 0)
  control[names(conlist)] <- conlist
  control <- list(maxit=as.integer(control$maxit),
                  tol=as.double(control$tol),
                  smooth=as.integer(control$smooth),
                  intsmooth=as.integer(control$smooth))
}

kern.str <- function(kernlist) {
  kernstr <- list(kern = 1, poly = 1, band = 1)
  kernstr[names(kernlist)] <- kernlist
  kernstr <- list(kern=as.integer(kernstr$kern),
                  poly=as.integer(kernstr$poly),
                  band=as.double(kernstr$band))
}
