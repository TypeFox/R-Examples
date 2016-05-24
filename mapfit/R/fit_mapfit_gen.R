
mapfit.gen.options <- function() {
  list(poisson.eps = sqrt(.Machine$double.eps),
       uniform.factor = 1.01,
       uniform.atol = sqrt(.Machine$double.eps),
       maxiter = 2000,
       reltol = sqrt(.Machine$double.eps),
       abstol = +Inf,
       annealing = FALSE,
       temperature = seq(0.9, 1, length.out=10),
       annealing.iter = NULL)
}

mapfit.gen.verbose <- function() {
  list(emstep = FALSE,
    emprogress = 1)
}

mapfit.gen <- function(map, data, initialize = TRUE, stationary = TRUE,
  control = list(), verbose = list(), ...) {
  call <- match.call()

  con <- mapfit.gen.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  ver <- mapfit.gen.verbose()
  nmsC <- names(ver)
  ver[(namc <- names(verbose))] <- verbose
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in verbose: ", paste(noNms, collapse = ", "))

  ## init parameters
  if (initialize) {
    map <- emfit.init(model=map, data=data, verbose=ver)
  }

  tres <- system.time(result <- emfit(map, data, initialize=FALSE,
    ufact=con$uniform.factor, eps=con$poisson.eps, atol=con$uniform.atol,
    control=con, verbose=ver, stationary=stationary, ...))
  result <- c(result, list(stationary=stationary, data=data@data, ctime=tres[1], call=call))
  class(result) <- "mapfit.result"
  result
}

print.mapfit.result <- function (x, ...) {
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
