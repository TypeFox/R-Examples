##
## File:   ProcessNetCDandH5FFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the NetCDF and HDF5 file processor.
##
processNetCDandH5FFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("nc", "h5"))) return(NULL)

  ## Create the results holder.
  results <- list()

  ## Attempt to open the file.
  header <- ncdf4::nc_open(path)

  ## Scan the dimensions
  for (index in 1:header$var) {

    ## Define the variable name.
    variable <- header$var[[index]]

    ## Setup the processing command.
    command <- paste("try (results$'", variable,
      "' <- ncdf4::ncvar_get(",
      "header, header$var[[index]]))", sep = "")

    ## Attempt to evaluate the processing command.
    eval(parse(text = command))

  }

  ## Return the results.
  results

}
