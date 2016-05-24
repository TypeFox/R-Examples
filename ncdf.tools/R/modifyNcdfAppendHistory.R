modifyNcdfAppendHistory = function(
    ##title<< Append a string to netCDF history
    file ##<< character sting or RNetCDF file connection: file to write to.
, string ##<< character string: string to append to the history
)
  ##description<<
  ## Convenience function to append a string together with the date and the user
  ## to the history attribute of an NetCDF file. 
  {
   if (class(file) == 'character') {
     con <- open.nc(file, write = TRUE)
   } else {
     con = file
   }
   history = ''
   if (is.element('history', infoNcdfAtts(con, 'NC_GLOBAL')[, 'name'])) {
     history <- att.get.nc(con, 'NC_GLOBAL', 'history')
   } 
   history.new <- paste(history,  Sys.time(), ':', string, 'by', Sys.info()['user'])
   att.put.nc(con, 'NC_GLOBAL', 'history', 'NC_CHAR', history.new)
   if (class(file) == 'character')
     close.nc(file)
   ##value<< nothing is returned
}
