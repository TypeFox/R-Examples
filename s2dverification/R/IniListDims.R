IniListDims <- function(dims, lenlist) {
  u <- list()
  for (jdim in 1:length(dims)) {
    u[[jdim]] <- 1:dims[jdim]
  }
  for (jdim in (length(dims) + 1):lenlist) {
    u[[jdim]] <- 1
  }
  #
  #  Outputs
  # ~~~~~~~~~
  #
  out <- u
}
