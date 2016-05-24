aggregateNcdf <- function(
  ##title<< Aggregate data in netCDF files
  fileName ##<< character vector: names of the files to aggregate.
  , path.out = getwd() ##<< character: path to save the results files to. 
  , period ##<< integer or one of hour, day, month or year: period to aggregate to. In case
           ##   of an integer value, the unit is time steps.
  )
  ##description<<
  ## This function aggregates time periods in netCDF files. Basically it is just a
  ## wrapper around the respective cdo function.
{

  ##test input
  if (!checkInstalled('cdo'))
    stop('cdo not found. Please install it.')

  if (inherits(period, 'character')) {
    cdoOperator <- paste(switch(period, hourly='hour', daily='day', monthly='mon',
                        yearly='year', 'not_available'), 'mean', sep = '')
    if(cdoOperator == 'not_vailable')
      stop(paste('Period \'', period, '\' not implemented!', sep = ''))
  } else if (inherits(period, 'numeric')) {
    cdoOperator      <- paste('timselmean,', period, sep = '') 
  } else {
    stop('Plase only supply numeric or character values for \'period\'!')
  }
  
  ## determine cdo command values
  fileNameOut <- file.path(path.out, sub('[.]nc', paste('_', period, '.nc', sep = ''), sub('.*/', '', fileName)))
  
  cdoCmd      <- paste('cdo ', cdoOperator, ' ', fileName, ' ', fileNameOut, sep = '')

  ##run aggregation
  system(cdoCmd)
  cat(paste('Created file ', fileNameOut, '.\n', sep = ''))
  ##value<<
  ## character string: name of the file created. 
  invisible(fileNameOut)
}


