

rvpar <- function (...) {
  args <- list(...)
  rvpar  <- getOption("rv")
  Pars <- rvpar$par
  if (length(args)==0) {
    return(Pars)
  }
  arg.names <- names(args)
  if (is.null(arg.names)) {
    args <- unlist(args)
    p <- Pars[args]
    if (length(args)==1) {
      return(p[[args]])
    } else {
      return(p)
    }
  }
  oldpar <- Pars
  for (name in arg.names) { 
    if (nchar(name)>=1) {
      Pars[name] <- args[name]
    }
  }
  rvpar$par <- Pars
  options(rv=rvpar)
  return(oldpar)
}


