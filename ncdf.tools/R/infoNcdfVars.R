infoNcdfVars <- function(
##title<< Display information about all variables in netCDF file
  file.con            ##<< a NetCDF object pointing to the respective netCDF file.
  , order.var = c('id', 'name')[2] ##<< character vector: 
                      ##   Whether to sort the variables according to their name (default) or id.
  , info.ext = FALSE  ##<< logical: whether to compute ranges/means etc. for the variables.
                      ##  Setting this to TRUE may take a while to compute with large files.
  , dimvars = FALSE   ##<< logical: whether to include the coordinate variables in the output. 
  
)
  ##description<< This function returns different summary information about all variables
  ##              in a netCDF file.
  ##seealso<<
  ##\code{\link{infoNcdfDims}}, \code{\link{infoNcdfAtts}}
{
  close.file =  FALSE
  if (inherits(file.con,  'character')) {
    if (!file.exists(file.con))
      stop(paste('File ', file.con, ' not existent!'))
      file.con <- open.nc(file.con)
    close.file = TRUE
  } else if (!class(file.con) == 'NetCDF') {
    stop(paste('file.con needs to be a connection to an open NetCDF file and of class',
               'NetCDF (library(RNetCDF))', sep=''))
  }
  
    #determine dim names and number of vars/dims
    n.vars    <- file.inq.nc(file.con)$nvars
    n.dims    <- file.inq.nc(file.con)$ndims
    dim.names <- character(length = n.dims)
    dim.ids   <-  0:(n.dims - 1)
    for (h in dim.ids)
        dim.names[h+1] <- dim.inq.nc(file.con, h)$name

    #determine characteristics for each single variable
    dim.names.var <- matrix(NA, n.vars, n.dims)
    colnames(dim.names.var)<- paste(1:n.dims, '.dim', sep = '')
    varname.var   <- character(length = n.vars)
    n.dims.var    <- integer(length = n.vars)
    units.var     <- character(length = n.vars)
    dimids.var    <- matrix(NA, n.vars, n.dims)
    colnames(dimids.var)<-paste('dim.id.', 1:n.dims, sep = '')
    type.var      <- character(length = n.vars)
    n.values.var  <- integer(length = n.vars)
    range.var     <- character(length = n.vars)

    #loop through all variables
    for (i in 1:n.vars) {
        varname.var[i]       <- var.inq.nc(file.con, i - 1)$name
        n.dims.var[i]        <- var.inq.nc(file.con, i - 1)$ndims
        type.var[i]          <- var.inq.nc(file.con, i - 1)$type
        dimids.var[i, 1:var.inq.nc(file.con, i - 1)$ndims] <- var.inq.nc(file.con, i - 1)$dimids
        if (info.ext) {
            data.dummy       <- var.get.nc(file.con, i - 1)
            range.var[i]     <- paste(round(range(data.dummy, na.rm = TRUE),
                                                digits = 2), collapse = '-')
            n.values.var[i]  <- sum(!is.na(data.dummy))
        }
        dims.cols            <- match(var.inq.nc(file.con, i - 1)$dimids, dim.ids)
        dim.names.var[i, dims.cols] <- dim.names[dims.cols]
        att.unit.id          <- grep('Unit', infoNcdfAtts(file.con, i - 1)[, 1], ignore.case = TRUE)
        if (length(att.unit.id) == 1)
            units.var[i]     <- att.get.nc(file.con, i-1, att.unit.id - 1)
    }

    #compile info into 1 dataframe
    ##value<<
    ## a dataframe with the different information in its columns and each variable in one row.
    varinfo.out <- data.frame(id=0:(n.vars - 1),
                              name = varname.var,
                              unit = units.var,
                              n.dims = n.dims.var,
                              type = type.var,
                              n.values = n.values.var,
                              range = range.var,
                              stringsAsFactors = FALSE)
    varinfo.out   <- cbind(varinfo.out, dim.names.var, dimids.var)
    col.order     <- match(order.var, colnames(varinfo.out))
    varinfo.out   <- varinfo.out[order(varinfo.out[, col.order]), ]
    rowsCoordVars <- na.omit(match(dim.names, varinfo.out$name))
    if (length(rowsCoordVars) > 0 & !dimvars)
      varinfo.out   <- varinfo.out[-rowsCoordVars,]
  if (close.file)
    close.nc(file.con)
  return(varinfo.out)
}




