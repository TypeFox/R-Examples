transNcdfSubset <- function(
  ##title<< Cut and save a subset of a netCDF file
  file.input ##<< character string: name of the input ncdf file. 
  , dim.values = list(latitude =c(), longitude = c(), time = c())
  , values.type =c('range', 'indices', 'values')[2] ##<< character string:
           ## type of the dim.values supplied. 'range' means that the lower an upper
           ## border are supplied, 'indices' means that 1:n indices are supplied,
           ## 'values' would imply actual coordinate values. 
  , file.output = sub('[.]nc', '_subs.nc', file.input) ##<< character string: name
           ## of the results file.
  , var.name = readNcdfVarName(file.input))
##description<<
## This function reads a subset of lat/lon/time values out of a netCDF file and creates
## a new netCDF file with the results.
{
  ##ToDo facilitate other scenarios than lat/lon/time
  ##TODO merge with transNcdfCutTimes

  #determine dimension indices
  dim.indices=list(latitude=c(), longitude=c(), timestep=c())
  con.source <- open.nc(file.input)
  dimvals.src        <- readNcdfCoordinates(con.source)
  if (sum(is.na(match(c('latitude', 'longitude', 'time'),  names(dimvals.src) ))) >  0)
    stop('Dimension names in file not latitude/longitude/time. Rename!')

  
  if (values.type == 'range') {
    for (dimT in names(dim.values))
      dim.indices[[dimT]] <-  which(dimvals.src[[dimT]] >= min(dim.values[[dimT]]) & dimvals.src[[dimT]] <= max(dim.values[[dimT]]))    
  } else  if (values.type == 'indices') {
    dim.indices = dim.values
  } else if (values.type == 'values') {
    for (dimT in names(dim.values)) 
      dim.indices[[dimT]] <- match(dim.values[[dimT]], dimvals.src[[dimT]])
  } else {
    stop(paste('Value for values.type of \'', values.type, '\' not supported!', sep = ''))
  }

  # create target file
  dimvals.trgt <-  list()
  for (dimT in names(dimvals.src))
    dimvals.trgt[[dimT]] <- dimvals.src[[dimT]][dim.indices[[dimT]]]
  con.target <- create.nc(file.output)
  for (dimT in names(dimvals.src)) {
    dim.def.nc(con.target,  dimT, dimlength = length(dimvals.trgt[[dimT]]))
    var.def.nc(con.target, dimT, 'NC_DOUBLE',  dimensions =  dimT)
    var.put.nc(con.target, dimT, dimvals.trgt[[dimT]])
    modifyNcdfCopyAtts(con.source, dimT, dimT, con.target)
  }
  modifyNcdfCopyAtts(con.source, 'NC_GLOBAL', 'NC_GLOBAL', con.target)

  # load, transpose and subset source data
  results.rotation  <- transNcdfRotate(con.source)
  data.orig <- results.rotation$data
  aperm.reverse <-  results.rotation$aperm.reverse
  dim.sub <- c(length( dim.values$latitude), length( dim.values$longitude),  length( dim.values$time))
  data.target.sub <- array(data.orig[indexDimVecs2Matrix(dim.indices$latitude, dim.indices$longitude, dim.indices$time)], dim = dim.sub)

  # write to target
  var.def.nc(con.target, var.name, var.inq.nc(con.source, var.name)$type, var.inq.nc(con.source, var.name)$dimids )
  modifyNcdfCopyAtts(con.source, var.name, var.name, con.target)
  
  var.put.nc(con.target, var.name, aperm(data.target.sub, perm = aperm.reverse))
  close.nc(con.target)
  close.nc(con.source)

  cat(paste('Finalized file ', file.output, '\n', sep = ''))
  ##value<< character string: name of the file created.
  invisible(file.output)
}


