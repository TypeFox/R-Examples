convertBinary2Ncdf = function(
##title<< transform binary file to netCDF file
##description<< This function transforms a binary data file to a netCDF file formatted
##              in a standardized way.
    file.input                             ##<< character string: name of the input file.
    , date.vec                             ##<< R date object: time vector for the time coordinate
    , length = 1                           ##<< integer: Length in bytes of each entry in the input file.
    , type = numeric()                     ##<< R data type of the data in the input file.
    , type.ncdf = 'NC_DOUBLE'              ##<< character string: Desired data type in the netCDF file.
    , dimensions                           ##<< character vector: Names of the dimensions in the binary file.
    , dimension.values                     ##<< list: Each list element has to contain the coordinate values for
                                           ##   the respective dimension.
    , signed = TRUE                        ##<< logical: Whether the binary file contains signed integer values.
    , var.name                             ##<< character string: Short name of the variable in the binary file
                                           ##   (used for the meta data in the NetCDF file).
    , long_name = var.name                 ##<< character string: long name of the variable in binary file
                                           ##   (used for the meta data in the NetCDF file).
    , var.units = '[]'                     ##<< character string: units of the variable
                                           ##   (used for the meta data in the NetCDF file).
    , scale.factor.in = 1                  ##<< numeric: factor to multiply the binary input data with.
    , scale.factor.out = scale.factor.in   ##<< numeric: desired scale factor of the data in the netCDF file.
    , na.value.in = -9999                  ##<< numeric: missing value for input data.
    , na.value.out = na.value.in           ##<< numeric: missing value for output data.
    , offset.in = 0                        ##<< numeric: offset for input data.
    , offset.out = offset.in)              ##<< numeric: offset for output data.
##value<< Nothing is returned but a netCDF file with a standardized name is written
##        in the working directory.
{
  call.args <- convertArgs2String()
  cat('Loading data ...\n')
  dims.lengths <- sapply(dimension.values,length)
  to.read      <- file(file.input, "rb")
  data.raw     <- readBin(to.read, type, n = prod(dims.lengths), endian = "little",
                          size=length, signed = signed)
  close(to.read)
  cat('Transforming data ...\n')
  data.array   <- array(data.raw,dim=dims.lengths)
  data.array[data.array==na.value.in]=NA
  data.array   <- (data.array * scale.factor.in/scale.factor.out) + (offset.in -offset.out)
  cat('Creating netCDF file ...\n')
  file.name    <- paste(var.name,dims.lengths[dimensions=='longitude'],
                        dims.lengths[dimensions=='latitude'],
                        as.character(years(min(date.vec))),
                        as.character(years(max(date.vec))),'nc',sep='.')
  createStdNcdfFile(file.name, var.name, lat.length = dims.lengths[dimensions=='latitude'],
                    long.length = dims.lengths[dimensions=='longitude'],
                    time.length = dims.lengths[dimensions=='time'],
                    add_offset = offset.out, scale_factor = scale.factor.out,
                    missing_value = na.value.out, type.var = 'NC_SHORT', units =  var.units)
  cat('Writing data ...\n')
  file.con      <- open.nc(file.name,write=TRUE)
  data.array.perm <- aperm(data.array,perm=match(infoNcdfDims(file.con)[,'name'],dimensions))
  var.put.nc(file.con,var.name,data.array)
  var.put.nc(file.con,'latitude',dimension.values[dimensions=='latitude'][[1]])
  var.put.nc(file.con,'longitude',dimension.values[dimensions=='longitude'][[1]])
  time.lilian <- as.numeric(julian(dimension.values[dimensions=='time'][[1]],
                                   origin = as.POSIXct("1800-01-01", tz="UTC")))
  var.put.nc(file.con,'time',time.lilian)
  modifyNcdfDefAtts(file.con,var.name,atts=list(long_name=long_name))
  history.string <- paste(Sys.time(),' netCDF file created from file ', file.input,
                          ' using function data.bin2ncdf() by ', Sys.info()['login'], sep = '')
  modifyNcdfDefAtts(file.con, 'NC_GLOBAL', atts = list(history = history.string,
                                             creation_settings = call.args))
  close.nc(file.con)
  cat('File sucessfully transformed!\n')
}
