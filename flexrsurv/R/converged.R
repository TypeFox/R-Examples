converged <- function(fit, step="init step", stop=FALSE){
  if (fit$convergence==1) {
    ret <- FALSE
    if (stop) stop(paste("iteration limit 'maxit' have been reached", step, sep=" "))
  }
  else if (fit$convergence==10) {
    ret <- FALSE
    if (stop) stop(paste("Degenerancy of Nelder_Mead simplex", step, sep=" "))
  } else {
    ret <- TRUE
  }
  ret
} 
