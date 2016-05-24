createLatLongTime <- function(
  ##title<< Create empty lat/lon/time netCDF file
  file.name                           ##<< character string: name of the target file.
  , var.names = sub('[.]nc', '', file.name) ##<< character vector: names of the variables
                                      ##         in the target file.
  , lat.values = c()                  ##<< numeric values: coordinate values for the latitude
                                      ##   positions.
  , long.values = c()                 ##<< numeric values: coordinate values for the latitude
                                      ##   positions.
  , time.values = c()                 ##<< POSIXct vector: time values for the time dimension
  , add.dims =  list()
  , lat.length  = length(lat.values)  ##<< integer: length of the latitude dimension
  , long.length = length(long.values) ##<< integer: length of the longitude dimension
  , time.length = length(time.values) ##<< integer: length of the time dimension
  , scale_factor = 1       ##<< numeric: scale factor
  , add_offset = 0         ##<< numeric: offset
  , type.var = 'NC_DOUBLE' ##<< character string: type of the data
  , missing_value = -9999  ##<< numeric: missing data value
  , units = '[]'           ##<< character string: units of the variables in target file.  
)
  ##description<<
  ## this function creates an empty standardized latitude/longitude/time netCDF file.
  ##value<<
  ## Nothing is returned but a file is created. 
{
  ##TODO: units has to work with more than one variable
  file.con  <- create.nc(file.name)
  dims.used = c()
  
  if (0 != lat.length) {
    dims.used = c(dims.used, 'latitude')
    dim.def.nc(file.con, 'latitude', dimlength = lat.length)
    var.def.nc(file.con, 'latitude','NC_DOUBLE', 'latitude')
    modifyNcdfDefAtts(file.con, 'latitude', atts = list(long_name = "latitude",
                                              units = "degrees_north" ,
            standard_name = "latitude"))       
    if (length(lat.values) > 0)
      var.put.nc(file.con, 'latitude', lat.values[order(lat.values, decreasing = TRUE)])  
  }   
  if (0 != long.length) {
    dims.used = c(dims.used, 'longitude')
    dim.def.nc(file.con, 'longitude', dimlength = long.length)
    var.def.nc(file.con, 'longitude','NC_DOUBLE', 'longitude')
    modifyNcdfDefAtts(file.con,'longitude',atts = list(long_name = "longitude",
                                             units = "degrees_east" ,
            standard_name = "longitude"))
    if (length(long.values) > 0)
      var.put.nc(file.con, 'longitude', long.values[order(long.values)])  
  }

  #define additional dimensions
  if (0 != length(add.dims)) {
    for (i in 1:length(add.dims)) {
      dims.use = names(add.dims)[i]
      start =  1
      count = length(add.dims[[i]])
      if (classR2Ncdf(add.dims[[i]]) == 'NC_CHAR') {
        if (!is.element('max_string_length', infoNcdfDims(file.con)$name))
          dim.def.nc(file.con, 'max_string_length', max(nchar(add.dims[[i]])))
        dims.use =  c('max_string_length', dims.use)
        start = c(1, start)
        count = c( max(nchar(add.dims[[i]])), count)
      }  
      dim.def.nc(file.con, names(add.dims)[i], dimlength = length(add.dims[[i]]))
      var.def.nc(file.con, names(add.dims)[i], classR2Ncdf(add.dims[[i]]),  dims.use)
      var.put.nc(file.con, names(add.dims)[i],  add.dims[[i]], start= start, count = count)
      dims.used <-  c(dims.used, names(add.dims)[i])
    }
  }

  if (0 != time.length) {
    dims.used =  c(dims.used, 'time')
    dim.def.nc(file.con, 'time', dimlength = time.length)
    var.def.nc(file.con, 'time','NC_DOUBLE', 'time')
    if (length(time.values) > 0) {
      if (inherits(time.values, c('POSIXlt', 'POSIXct'))) {
        modifyNcdfDefAtts(file.con,'time',atts = list(long_name = "time",units = "days since 1800-01-01 00:00" ,
                                            calendar = "gregorian"))     
        time.values <- as.numeric(julian(time.values, origin = as.POSIXct("1850-01-01", tz="UTC")))
      }
      var.put.nc(file.con, 'time', time.values)
    }     
  }
  
  # define attributes
  for (var.name.t in var.names) {
    var.def.nc <- var.def.nc(file.con, var.name.t, type.var, dims.used)
    modifyNcdfDefAtts(file.con, var.name.t, atts = list(scale_factor = scale_factor,
                                              add_offset = add_offset,
                                              missing_value = missing_value,
                                              `_FillValue` = missing_value, units = units))
  }
  hist_string <- paste('File created on ', Sys.time(), ' by ', Sys.info()['user'] , sep = '')
  att.put.nc(file.con, 'NC_GLOBAL', 'history', 'NC_CHAR', hist_string)
  close.nc(file.con)
  cat(paste('Created file', file.name), '\n')
}
