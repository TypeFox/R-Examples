modifyNcdfAddDim <- function(
##title<< Add a new dimension to one or more variables in a netCDF file
        file.con.orig          ##<< a NetCDF object pointing to the respective netCDF file FROM which to copy
        , file.con.copy        ##<< a NetCDF object pointing to the respective netCDF file TO which to copy
        , var.name = 'Default' ##<< character vector: names of the variables to which a dimension should be added. Defaults to
                               ##   all except those with identical names as dimensions in file.con.orig (coordinate variables)
        , dim.name = 'new.dim' ##<< character string: name of the dimension to add
        , dim.values = c()     ##<< numeric/character vector with the values for the dimension (coordinate values)
        , dim.length = length(dim.values)  ##<< integer: length of the dimension to add
        , dim.pos.copy = 1     ##<< integer: position in the new dimension where to copy the original data. If set to 0,
                               ##   no values are copied and the variable in the new file will be empty. Setting to
                               ##   1 (default) results in the original values to be filled in the first value of
                               ##   the new dimension and the remaining values left empty (NaN).
)
  ##description<<
  ## Adds another dimension to specified variables in a netCDF file and saves the results in
  ## another netCDF file.
  ##seealso<< 
  ##\code{\link{modifyNcdfCopyMetadata}}, \code{\link[RNetCDF]{att.copy.nc}},
  ## \code{\link{modifyNcdfCopyVar}}
{
    #test input
    if (var.name == 'Default') {
        var.name = setdiff(infoNcdfVars(file.con.orig)$name, infoNcdfDims(file.con.orig)$name)
    } else if (sum(!is.element(var.name, infoNcdfVars(file.con.orig)$name)) > 0) {
        stop('Variable name(s) not available in file.con.orig!')
    }
    if(dim.length == 0)
        stop('Dimension of length Zero makes no sense!')
    if((!length(dim.values) == 0 ) && (!dim.length == length(dim.values)))
        stop('Wrong dim.length specified!')

    #copy attributes
    modifyNcdfCopyMetadata(file.con.orig, file.con.copy)

    #define new dimension
    dim.def.nc(file.con.copy, dim.name, dim.length, unlim=FALSE)
    if (length(dim.values) > 0) {
        type.dim<- classR2Ncdf(dim.values)
        var.def.nc(file.con.copy, dim.name, type.dim, dim.name)
        var.put.nc(file.con.copy, dim.name, dim.values)
    }
    id.new.dim  <- dim.inq.nc(file.con.copy, dim.name)$id
    
    #create variable
    var.name    <- sort(var.name)
    vars.orig   <- infoNcdfVars(file.con.orig)[match(var.name, infoNcdfVars(file.con.orig)$name), ]
    for (i in 1:length(var.name)) {
        print(paste('Add netCDF dimension: Processing variable ', i, ' of ', length(var.name), ' : ', var.name[i], sep=''))
        vars.orig.dims.t <- var.inq.nc(file.con.orig, var.name[i])$dimids
        vars.copy.dims.t <- c(vars.orig.dims.t, id.new.dim)
        var.def.nc(file.con.copy, vars.orig[i, 'name'], vars.orig[i,'type'], vars.copy.dims.t)
        modifyNcdfCopyAtts(file.con.orig, var.copy = var.name[i], var.name[i], file.con.copy)
        #copy values
        if (!(dim.pos.copy == 0)) {
            data.transfer    <- var.get.nc(file.con.orig, var.name[1])
            start.ind <- rep(1, times = length(vars.copy.dims.t))
            stop.ind  <- infoNcdfDims(file.con.copy)$length
            start.ind[id.new.dim + 1] <-  dim.pos.copy
            stop.ind[id.new.dim + 1]  <-  dim.pos.copy
            var.put.nc(file.con.copy, var.name[i], data.transfer, start = start.ind, count = stop.ind)
        }
        sync.nc(file.con.copy)
    }

}
