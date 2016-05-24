evmSimSetSeed <- function(x){
  if (oldClass(x) != "evmSim"){
      stop("This function expects an object of class \'evmSim\'.")
  }

  assign(".Random.seed", x$seed, envir=.GlobalEnv)

  invisible(x$seed)
}
