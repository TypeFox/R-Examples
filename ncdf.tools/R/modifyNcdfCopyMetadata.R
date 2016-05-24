modifyNcdfCopyMetadata <- function(
##title<< Copy attributes and dimensions between netCDF files
    file.con.orig        ##<< a NetCDF object pointing to the respective netCDF file from which to copy
    , file.con.copy      ##<< a NetCDF object pointing to the respective netCDF file to which to copy
    , glob.atts = TRUE   ##<< logical: whether to copy all global attributes
    , dimensions = TRUE  ##<< logical: whether to copy all dimensions
    )
  ##description<<
  ## This function copies all global attributes and/or all dimensions from one
  ## netCDF file to another.

  ##seealso<<
  ##\code{\link{modifyNcdfCopyVar}}, \code{\link[RNetCDF]{att.copy.nc}}

{
  n.dims <- file.inq.nc(file.con.orig)$ndims
  n.atts <- file.inq.nc(file.con.orig)$ngatts
  
  if (glob.atts) {
    for (i in 1:n.atts){
      att.name  <- att.inq.nc(file.con.orig, 'NC_GLOBAL', i - 1)$name
      att.type  <- att.inq.nc(file.con.orig, 'NC_GLOBAL', i - 1)$type
      att.value <- att.get.nc(file.con.orig, 'NC_GLOBAL', att.name)
      att.put.nc(file.con.copy, 'NC_GLOBAL', att.name, att.type, att.value)
    }
  }
  if (dimensions) {
    for (i in 1:n.dims) {
      dim.name   <- dim.inq.nc(file.con.orig, i - 1)$name
      dim.length <- dim.inq.nc(file.con.orig, i - 1)$length
      dim.def.nc(file.con.copy, dim.name, dim.length)
      if (is.element(dim.name, infoNcdfVars(file.con.orig, dimvars = TRUE)$name)) {
        var.def.nc(file.con.copy, dim.name, 'NC_DOUBLE', i - 1)
        var.put.nc(file.con.copy, dim.name, var.get.nc(file.con.orig, dim.name))  
        modifyNcdfCopyAtts(file.con.orig, dim.name, dim.name, file.con.copy)
      }
    }
  }
  ##value<< nothing is returned.
}
