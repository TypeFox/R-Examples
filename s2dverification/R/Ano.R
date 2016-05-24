Ano <- function(var, clim) {
  #
  #  Duplicate clim dimensions to heve same dimensions as var
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimvar <- dim(var)
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[2] != dimvar[2])) {
    clim <- InsertDim(clim, 2, dimvar[2])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[3] != dimvar[3])) {
    clim <- InsertDim(clim, 3, dimvar[3])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[4] != dimvar[4])) {
    clim <- InsertDim(clim, 4, dimvar[4])
  }

  #
  #  Raw anomalies
  # ~~~~~~~~~~~~~~~
  #
  ano <- var - clim

  #
  #  Outputs
  # ~~~~~~~~~
  #
  ano
}
