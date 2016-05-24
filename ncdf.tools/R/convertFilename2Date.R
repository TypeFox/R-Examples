convertFilename2Date <- function(
  ##title<< Convert file name patterns to R date object
    file.names        ##<< character vector: names of the files
    , fun.extr.string ##<< function
    , fun.conv.string ##<< function
  )
##description<<
## This function converts parts of netCDF file names to date strings (in case the
## file name contains date information). This is used, e.g.
## in transNcdfCutFiles.   
##seealso<<
##\code{\link{transNcdfCutFiles}}  
{
   date   <- do.call(fun.conv.string, list( do.call(fun.extr.string, list(file.names))))
   ##value<< POSIXct object with the date.
   return(date)
}
