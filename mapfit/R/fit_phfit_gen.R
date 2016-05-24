######## PHFIT for general PH

phfit.gen.options <- function() {
  list(poisson.eps = sqrt(.Machine$double.eps),
       uniform.factor = 1.01,
       maxiter = 2000,
       maxiter.init = 5,
       reltol = sqrt(.Machine$double.eps),
       abstol = +Inf,
       diff.init = c(1, 4, 16, 64, 256, 1024),
       scale.init = c(0.5, 1.0, 2.0),
       annealing = FALSE,
       temperature = seq(0.9, 1, length.out=10),
       annealing.iter = NULL)
}

phfit.gen.verbose <- function() {
  list(emstep = FALSE,
    emprogress = 1,
    cf1init = TRUE)
}

phfit.gen <- function(ph, data, initialize = TRUE, control = list(), verbose = list(), ...) {
  call <- match.call()

  con <- phfit.gen.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  # if (con$annealing.iter == NULL) {
  #   con$annealing.iter <- rep(5, length(con$temperature))
  # }

  ver <- phfit.gen.verbose()
  nmsC <- names(ver)
  ver[(namc <- names(verbose))] <- verbose
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in verbose: ", paste(noNms, collapse = ", "))

  ## init parameters
  if (initialize) {
    ph <- emfit.init(model=ph, data=data, verbose=ver)
  }

  tres <- system.time(result <- emfit(ph, data, initialize=FALSE,
    ufact=con$uniform.factor, eps=con$poisson.eps, control=con, verbose=ver, ...))
  result <- c(result, list(data=data@data, ctime=tres[1], call=call))
  class(result) <- "phfit.result"
  result
}

print.phfit.result <- function (x, ...) {
  cat("\n")
  cat(sprintf("Maximum LLF: %f\n", x$llf))
  cat(sprintf("AIC: %f\n", x$aic))
  cat(sprintf("Iteration:  %d / %d\n", x$iter, x$control$maxiter))
  cat(sprintf("Computation time (user): %f\n", x$ctime))
  cat(sprintf("Convergence: %s\n", x$convergence))
  cat(sprintf("Error (abs): %e (tolerance %e)\n", x$aerror, x$control$abstol))
  cat(sprintf("Error (rel): %e (tolerance %e)\n", x$rerror, x$control$reltol))
  cat("\n")
  emfit.print(x$model)
  cat("\n\n")
  invisible(x)
}
