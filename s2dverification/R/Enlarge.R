Enlarge <- function(var, numdims) {
  #
  #  Enlarge the number of dimensions to 20 --> enlvar
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  d <- c(dimsvar, array(1, dim = 20))
  enlvar <- array(dim = d[1:20])
  enlvar[, , , , , , , , , , , , , , , , , , , ] <- var
  #
  #  Reduce the number of dimensions to the required one
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  outvar <- array(dim = d[1:numdims])
  outvar[] <- enlvar
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar
}
