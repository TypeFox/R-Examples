readNcdfVarName <- function(
  ##title<< Get name of variable in netCDF file
  file ##<< connection to the netCDF file.
)


##description<<
## readNcdfVarName tries to automatically detect the name of the "main" variable in a netCDF file. The name returned is the
## name of a non coordinate variable. If more than one of such variables are existent, the name of the variable
## which spans all available dimensions or with a name appearing as a pattern in the file name is used.

##seealso<<
## \code{\link[RNetCDF]{RNetCDF}}, \code{\link{infoNcdfVars}}

{
  if (class(file) == 'character') {
    file.name <- file
    file.con <- open.nc(file)
  } else {
    file.con <- file
    file.name = ''
  }  
  var.name         <- infoNcdfVars(file.con, order.var ='id',  dimvars = FALSE)$name

  # exclude reserved names
  names.excluded   <- c("lon_bnds", "lat_bnds", "time_bnds", "borders.low", "borders.up" )
  var.name         <- setdiff(var.name, names.excluded)
  var.name         <- var.name[!grepl('flag.orig$', var.name)]

  #try to determine variable name based on matches in file name
  if (length(var.name) > 1 && nchar(file.name) > 0) {
    names.matched <- pmatch(var.name,  file.name)
    if (sum(!is.na(names.matched)) == 1) {
      var.name <-  var.name[!is.na(names.matched)]
    }
  }

  #try to determine variable name based dimension size matches  
  if(length(var.name) > 1) {
    vars.ndims <-  infoNcdfVars(file.con)$n.dims
    n.dims <- dim(infoNcdfDims(file.con))[1]
    if (sum(vars.ndims == n.dims ) == 1) {
      var.name <-  var.name[vars.ndims == n.dims ]
    }
  }
     
  if ((length(var.name) > 1) ) {
    stop('Not possible to detect variable name!')
  } 
  
  if (class(file) == 'character') {
    close.nc(file.con)
  } 
  ##value<< character string: name of the variable.   
  return(var.name)
}  
