
mapfit.gmmpp <- function(map, data, initialize = TRUE, stationary = TRUE,
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
