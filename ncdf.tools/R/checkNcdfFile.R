checkNcdfFile <- function(
##title<< check netCDF file for consistency with CF/COARDS/BGI netCDF conventions
  file.name     ##<<character string: file name to check
  , dims = c('longitude', 'latitude', 'time') ##<< vector of strings:
                ## names of the dimensions which need to be in the file.
  , type = 'strict'  ##<< character string:
                ##   if 'strict', then all aspects are checked. If this
                ##   is any other value, only aspects relevant for the processing of
                ##   decomp.ncdf are checked.
  , var.check ='single' ##<< character string: If 'single', then readNcdfVarName
                ## is used to infer the name of the variable in the target file
                ## which will then be checked,
)
  ##description<< This function checks whether a netCDF file is consistent with the parts of the COARDS/CF
  ##              netCDF conventions used in the BGI department (MPI for Biogeochemistry,  Jena,  Germany).
{
  con.check     <- open.nc(file.name)

  dims.nonvalid <- is.na(match(dims, c('longitude', 'latitude', 'time')))
  if (sum(dims.nonvalid) > 0)
    stop(paste('Function is not designed to check dimensions ', dims[dims.nonvalid], '!', sep=''))

  #check dimensions
  dims.exists <- match(dims, infoNcdfDims(con.check, extended = FALSE)[, 'name'])
  if (sum(is.na(dims.exists))>0) {
    close.nc(con.check)
    cat(paste('Dimension ', dims[is.na(dims.exists)], ' not existent!\n', sep=''))
    return(invisible(FALSE))
  }
  if (type == 'strict') {
    #check coordinate variables
    dims.exists <- match(dims, infoNcdfVars(con.check)[, 'name'])
    if (sum(is.na(dims.exists)) > 0)  {
      close.nc(con.check)
      cat(paste('Coordinate variable for dim ', dims[is.na(dims.exists)], ' not existent!\n', sep = ''))
      return(invisible(FALSE))
    }
    for (dim.t in dims) {
      if (!(sum(!is.na(var.get.nc(con.check, dim.t))) == dim.inq.nc(con.check, dim.t)$length)) {
        close.nc(con.check)
        cat(paste('Coordinate variable for dim ', dim.t, ' has missing values!\n', sep = ''))
        return(invisible(FALSE))
      }
    }
  } 

  #check file name
  if (type == 'strict') {
    var.name <- sub('[.].*', '', file.name)
    if (!is.element(var.name, infoNcdfVars(con.check)[,'name']))
      stop(paste('File does not have the main variable which has to be named ', var.name,
                 '(according to the file name targetvar.NumberLatitudes.NumberLongitudes.Year.nc )', sep = ''))    
  }
  
  #check attributes
  variables <- infoNcdfVars(con.check)[, 'name'][is.na(match(infoNcdfVars(con.check)[, 'name'], infoNcdfDims(con.check, extended = FALSE)[, 'name']))]
  variables <- variables[!is.element(variables, c("flag.orig", 'borders.low', 'borders.up'))]
  if (length(variables) > 1 && var.check == 'single')
    variables  <- readNcdfVarName(file.name)
  
  if (type == 'strict') {       
    atts.check <- c('_FillValue', 'missing_value', 'add_offset', 'scale_factor', 'units')
    for (var.t in variables) {
      atts.found <- match(atts.check, infoNcdfAtts(con.check, var.t)[, 'name'])
      if (sum(is.na(atts.found)) > 0) {
        close.nc(con.check)
        cat(paste('Variable ', var.t, ' does not have attributes ', paste(atts.check[is.na(atts.found)], collapse=',  '), '\n', sep=''))
        return(invisible(FALSE))
      }
    }
  }
  
  #check time vector
  if (type == 'strict') {
    if (is.element('time', dims)) {
      if (!is.element('units', infoNcdfAtts(con.check, 'time')[, 'name'])) {
        close.nc(con.check)
        cat('Time variable needs units attribute!\n')
        return(invisible(FALSE))
      }
      att.time.units <- att.get.nc(con.check, 'time', 'units')
      if (!(att.time.units == 'days since 1800-01-01 00:00')) {
        close.nc(con.check)
        cat('Change time vector to days since 1800-01-01 00:00 !\n')
        return(invisible(FALSE))
      }
    }
  }
  #check missing_values attribute
  atts.check <- c('_FillValue', 'missing_value')
  for (var.t in variables) {
    if (sum(is.na(match(atts.check, infoNcdfAtts(con.check, var.t)[, 'name']))) == 2) {
      cat(paste('One of the attributes _FillValue or missing_value has to be available for variable ', var.t, '!.\n', sep=''))
      return(invisible(FALSE))
    }
    for (att.t in atts.check) {
      if (!is.na(match(att.t, infoNcdfAtts(con.check, var.t)[, 'name']))) {
        if (!(att.inq.nc(con.check, var.t, att.t)$type == var.inq.nc(con.check, var.t)$type)) {
          close.nc(con.check)
          cat(paste('Attribute ', att.t, ' of variable ', var.t, ' needs to be of the same class as ', var.t, '.\n', sep=''))
          return(invisible(FALSE))
        }
      }
    }
  }
  #check non obligatory global atts
  if (type == 'strict') {
    atts.check.global <- c('title', 'reference', 'history', 'provided_by', 'created_by')
    atts.missing      <- is.na(match(atts.check.global, infoNcdfAtts(con.check, 'NC_GLOBAL')[, 'name']))
    if (sum(atts.missing)>0)
      cat(paste('Consider adding the following global attributes: ', paste(atts.check.global[atts.missing], collapse=',  '), '\n', sep=''))
  }
  close.nc(con.check)
  cat('File check passed!\n')
  ##value<< logical: (invisible) whether the file passed the check
  return(invisible(TRUE))
}
