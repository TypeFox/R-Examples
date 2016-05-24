readNcdfSpectral <- function(
  ##title<< Read the results of a spectral decomposition (from a netCDF file)
  fileName       ##<< character string: name of the netCDF file
  , varName      ##<< character string: name of the variable to extract.
  , rangeBandsGet##<< vector: Vector defining the bands to extract. Can be either
                 ##   logical with one TRUE/FALSE per band in the file or a numeric
                 ##   vector of length two with the lower and the upper spectral
                 ##   border.
  )
  ##description<<
  ## readNcdfSpectral reads spectrally decomposed ncdf data (i.e. the output of a call to decomposeNcdf).
{
  .funSum    <- function(x) {
    xUse <- x[bandsTake]
    if (sum(!is.na(xUse)) == 0) {
      NA
    } else {
      sum(xUse, na.rm = TRUE)
    }
  }
  conT <- open.nc(fileName)
  dataIn <- var.get.nc(conT, varName)
  bordersLow <- var.get.nc(conT, 'borders.low')
  bordersUp  <- var.get.nc(conT, 'borders.up')
  close.nc(conT)
  if (inherits(rangeBandsGet, 'logical')) {
    if (length(rangeBandsGet != length(bordersUp)))
      stop('rangeBandsGet needs to be of the same length as dimension spectral bands in case of logicals!')
    bandsTake = rangeBandsGet
  } else {
    bandsTake   <- bordersUp <= rangeBandsGet[2] & bordersLow >= rangeBandsGet[1]
  }
  if (sum(bandsTake) == 0) {
    stop('no valid band selected!')
  } else if (sum(bandsTake) == 1){
    dataSummed <- dataIn[,bandsTake]
  } else if (sum(bandsTake) > 1) {
    dimsSum    <- 1:(length(dim(dataIn))-1)
    dataSummed <- apply(dataIn, dimsSum, .funSum)
  }
  ##value<<
  ## matrix: the spectral bands defined. 
  return(dataSummed)
}
