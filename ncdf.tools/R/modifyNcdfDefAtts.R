modifyNcdfDefAtts <- structure(function(
##title<< Define a set netCDF attributes at once
    file.con      ##<< a NetCDF object pointing to the respective netCDF file.
    , var.id      ##<< the variable id (integer) or name (string) for which to define attributes.
    , atts        ##<< list: the attributes to define (see details or an example).
    )
##description<<
## Easily define a couple of attributes for a single netCDF variable in one step.
##details<<
##The atts attribute should be a list with as many elements as attributes should be added to the
##variable in the netCDF file. The names of the attributes are taken from the names of the
##elements of this list and the attribute values are defined by the values of the list elements.
##The type/class of the attribute (values) is determined automatically.

##seealso<<
##\code{\link[RNetCDF]{att.put.nc}}

{
  n.steps              <- length(atts)
  att.names            <- names(atts)
  ind.fillvalue        <- na.omit(match(c('_FillValue', 'missing_value'), att.names))
  if (length(ind.fillvalue)==1) {
    att.t            <- c('_FillValue', 'missing_value')[!is.element(c('_FillValue', 'missing_value'), att.names)]
    att.names        <- c(att.names, att.t)
    atts             <- c(atts, atts[[ind.fillvalue]])
    names(atts)[length(atts)] <- att.t
    n.steps          <- n.steps + 1
  }
  for (i in 1:n.steps) {
    att.value        <- atts[[i]]
    if (is.numeric(att.value)) {
      if (att.names[i] == '_FillValue' | att.names[i] == 'missing_value') {
        att.type <- var.inq.nc(file.con, var.id)$type
      } else {
        att.type <- 'NC_DOUBLE'
      }
    } else {
      att.type     <- 'NC_CHAR'
    }
    att.put.nc(file.con, var.id, att.names[i], att.type, att.value)
  }
}, ex = function() {
  ## needs an open connection to a valid netCDF file pointed to by file.con
  attributes.define <- list(LongName = 'This is the long name',
                            missingValue = -99999,
                            units = 'm/s')
  library(RNetCDF)
  file.con   <- create.nc('test.nc')
  dim.def.nc(file.con, 'testdim')
  var.def.nc(file.con, 'test', 'NC_CHAR', 'testdim')
  modifyNcdfDefAtts(file.con, 'test', atts = attributes.define)

  ## show all attributes
  infoNcdfAtts(file.con, 'test')
})
