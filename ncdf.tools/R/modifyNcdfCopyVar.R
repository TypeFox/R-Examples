modifyNcdfCopyVar <- function(
##title<< Copy variable values between netCDF files
##description<< modifyNcdfCopyVar copies all values of one variable from one netCDF file to another netCDF file and
## takes care of dimensions etc. .
    file.con.orig        ##<< a NetCDF object pointing to the original netCDF file FROM which to copy the variable.
    , file.con.copy = file.con.orig     ##<< a NetCDF object pointing to the netCDF file TO which to copy the variable.
    , var.id.orig                      ##<< character string or netCDF variable id: The name or id of the variable to copy from.
    , var.id.copy = var.id.orig        ##<< character string or netCDF variable id: The name or id of the variable to copy to.
    )
##details<<
## Two cases are implemented: 
##
##Case 1: copy of one variable and attributes from one file to another file:
##The dimensions of the variable to copy have to be also existent (i.e. dimensions with the
##same name (not necessarily id)) in the netCDF file to which the variable
##should be copied. In addition these dimensions have to have the same sizes.
##
##Case 2: copy of one variable to another one (of different name) in the same file.
##seealso<<
##\code{\link{modifyNcdfCopyMetadata}}, \code{\link[RNetCDF]{att.copy.nc}}
{
    ##get values to copy from
    var.name.orig      <- var.inq.nc(file.con.orig, var.id.orig)$name
    var.type.orig      <- var.inq.nc(file.con.orig, var.id.orig)$type
    dim.names.file.orig<- infoNcdfDims(file.con.orig)$name
    dim.ids.var.orig   <- var.inq.nc(file.con.orig, var.id.orig)$dimids
    dim.lengths.orig   <- infoNcdfDims(file.con.orig)$length[dim.ids.var.orig + 1]
    
    ##get values to copy to
    if ((var.id.orig == var.id.copy) && (file.con.copy != file.con.orig)) {
      var.name.copy       <- var.name.orig
      dim.names.file.copy <- infoNcdfDims(file.con.copy)$name
      dim.ids.var.copy    <- match(dim.names.file.orig, dim.names.file.copy)-1
    } else if ((var.id.orig != var.id.copy) && (file.con.copy == file.con.orig)) {
      var.name.copy    <- var.id.copy
      dim.ids.var.copy <- dim.ids.var.orig
    } else {
      stop('Either var.ids or file connections (or both) have to be different!')
    }
    
    ##write to file
    var.def.nc(file.con.copy, var.name.copy, var.type.orig, dim.ids.var.copy)
    modifyNcdfCopyAtts(file.con.orig, var.name.orig, var.name.copy, file.con.copy)
    ##value<< 
    ##Nothing is returned
}
