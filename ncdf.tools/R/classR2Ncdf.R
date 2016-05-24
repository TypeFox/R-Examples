classR2Ncdf = function(
 object ##<< object which class should be determined
)
##title<<
## transfers R classes to netCDF classes

##description<<
## Crudely determines the netCDF class from R classes. Only integer, character
## and double are implemented yet.

##seealso<<
## \code{\link[RNetCDF]{RNetCDF}}
{
    if (class(object) == 'numeric' || class(object) == 'integer') {
        if ((sum(abs(object - round(object)) < .Machine$double.eps ^ 0.5) == length(object)) &
            (min(object, na.rm = TRUE) > -32.768 & max(object, na.rm = TRUE) < 32.767)) {
            ncdf.class <- 'NC_INT'
        } else {
            ncdf.class <- 'NC_DOUBLE'
        }
    } else {
        ncdf.class     <- 'NC_CHAR'
    }
    ##value<< character string: netCDF class used in the RNetCDF package.
    return(ncdf.class)
}
