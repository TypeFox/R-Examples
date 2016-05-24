#' Get variable name of the NetCDF file.
#' 
#' Get variable name in the NetCDF file. After knowning the name, you can use \code{loadNcdf} to load
#' the target variable.
#' 
#' @param filePath A path pointing to the netCDF file.
#' @return The names of the varialbes in the file.
#' @examples 
#' # First open the test NETcDF file.
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")
#' 
#' # Then if you don't know the variable name, you can use \code{getNcdfVar} to get variable name
#' varname <- getNcdfVar(filePath)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @import ncdf4
#' @references 
#' 
#' \itemize{
#' \item David Pierce (2015). ncdf4: Interface to Unidata netCDF (Version 4 or
#' Earlier) Format Data Files. R package version 1.14.1.
#' http://CRAN.R-project.org/package=ncdf4
#' }
#' 
#' 
#' 
#' @export
getNcdfVar <- function(filePath) {
  nc <- nc_open(filePath)
  names <- names(nc$var)
  return(names)
}


#' Load NetCDF file
#' 
#' @param filePath A path pointing to the NetCDF file, version3.
#' @param varname A character representing the variable name, you can use \code{getNcdfVar} to
#' get the basic information about the variables and select the target.
#' @param tz A string representing the time zone, default is GMT, if you know what time zone is 
#' you can assign it in the argument. If \code{tz = ''}, current time zone will be taken.
# @param drop When the time dimension only have one value, the output data will drop
# this dimension automatically (\code{drop = TRUE}), default value is \code{drop = FALSE}, then time dimension will be added.
# This argument mainly applies to the later calculations based on hyfo file. If the dimension
# is dropped, than some calculations may not be processed afterwards. 
#' @param ... Several arguments including Year, month, lon, lat 
#' type in \code{?downscaleNcdf} for details.You can load while downscale, 
#' and also first load than use \code{downscaleNcdf} to downscale.
#' @return A list object from \code{hyfo} containing the information to be used in the analysis, 
#' or biascorrection.
#' @examples 
#' # First open the test NETcDF file.
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")
#' 
#' # Then if you don't know the variable name, you can use \code{getNcdfVar} to get variable name
#' varname <- getNcdfVar(filePath)
#' 
#' nc <- loadNcdf(filePath, varname)
#' 
#' # you can directly add your downscale information to the argument.
#' nc1 <- loadNcdf(filePath, varname, year = 2006, lon = c(-2, -0.5), lat = c(43.2, 43.7))
#' nc2 <- loadNcdf(filePath, varname, year = 2005, month = 3:8, lon = c(-2, -0.5), 
#' lat = c(43.2, 43.7))
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @import ncdf4
#' @references 
#' 
#' \itemize{
#' \item David Pierce (2015). ncdf4: Interface to Unidata netCDF (Version 4 or
#' Earlier) Format Data Files. R package version 1.14.1.
#' http://CRAN.R-project.org/package=ncdf4
#' 
#' \item Santander MetGroup (2015). ecomsUDG.Raccess: R interface to the ECOMS User Data Gateway. R package
#' version 2.2-6. http://meteo.unican.es/ecoms-udg
#' }
#' 
#' 
loadNcdf <- function(filePath, varname, tz = 'GMT', ...) {
  nc <- nc_open(filePath)
  
  var <- nc$var
  # Use name to locate the variable
  call_1 <- as.call(c(
    list(as.name('$'), var, varname)
  ))
  var <- eval(call_1)
  if(is.null(var)) stop('No such variable name, check source file.')
  
  dimNames <- unlist(lapply(1:length(var$dim), function(x) var$dim[[x]]$name))
  
  # Only deals with the most common dimensions, futher dimensions will be added in future.
  dimIndex <- grepAndMatch(c('lon', 'lat', 'time', 'member'), dimNames)
  if (length(dimIndex) < 3) stop('Your file has less than 3 dimensions.')
  
  # First needs to identify the variable name, load the right data
  message('Loading data...')
  nc_data <- ncvar_get(nc, var)
  message('Processing...')
  
  gridData <- list()
  gridData$Variable$varName <- varname
  gridData$xyCoords$x <- var$dim[[dimIndex[1]]]$vals
  attributes(gridData$xyCoords$x)$name <- dimNames[dimIndex[1]]
  
  gridData$xyCoords$y <- var$dim[[dimIndex[2]]]$vals
  attributes(gridData$xyCoords$y)$name <- dimNames[dimIndex[2]]
  
  # Time part needs to be taken seperately
  
  timeUnit <- strsplit(var$dim[[dimIndex[3]]]$units, split = ' since')[[1]][1]
  timeDiff <- var$dim[[dimIndex[3]]]$vals
  # To get real time, time since when has to be grabbed from the dataset.
  timeSince <- as.POSIXlt(strsplit(var$dim[[dimIndex[3]]]$units, split = 'since')[[1]][2], tz = tz)
  
  
#  Date <- rep(timeSince, length(timeDiff))
  
  
  unitDic <- data.frame(weeks = 'weeks', days = 'days', hours = 'hours',
                        minutes = 'mins', seconds = 'secs')
  
  timeDiff <- as.difftime(timeDiff, units = as.character(unitDic[1, timeUnit]))
  
#   if (grepl('day', timeUnit)) {
#     Date$mday <- Date$mday + timeDiff
#   } else if (grepl('second', timeUnit)) {
#     Date$sec <- Date$sec + timeDiff
#   }
  Date <- timeSince + timeDiff
  
  # data directly loaded from ncdf4 will drop the dimension with only one value.
  # the varsize shows the real dimension, without any dropping.
  dim(nc_data) <- var$varsize 
  
  # Right now there is no need to add end Date, in furture, may be added as needed.
  gridData$Dates$start <- as.character(Date)
  
  # Assing data to grid data
  # At leaset should be 3 dimensions, lon, lat, time. So if less than 3, 
  
  gridData$Data <- nc_data
  
  attributes(gridData$Data)$dimensions <- dimNames
  
  if (!is.na(dimIndex[4])) gridData$Members <- var$dim[[dimIndex[4]]]$vals
  
  gridData$Loaded <- 'by hyfo package, http://yuanchao-xu.github.io/hyfo/'
  nc_close(nc)
  
  output <- downscaleNcdf(gridData, ...)
  
  return(output)
  
}




#' Downscale NetCDF file
#' @param gridData A hyfo list file from \code{\link{loadNcdf}}
#' @param year A vector of the target year. e.g. \code{year = 2000}, \code{year = 1980:2000}
#' @param month A vector of the target month. e.g. \code{month = 2}, \code{month = 3:12}
#' @param lon A vector of the range of the downscaled longitude, should contain a max value
#' and a min value. e.g. \code{lon = c(-1.5, 2,5)}
#' @param lat A vector of the range of the downscaled latitude, should contain a max value
#' and a min value. e.g. \code{lat = c(32,2, 36)}
#' @return A downscaled hyfo list file.
#' @examples 
#' # First open the test NETcDF file.
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")
#' 
#' 
#' # Then if you don't know the variable name, you can use \code{getNcdfVar} to get variable name
#' varname <- getNcdfVar(filePath)
#' 
#' nc <- loadNcdf(filePath, varname)
#' 
#' # Then write to your work directory
#' 
#' nc1 <- downscaleNcdf(nc, year = 2006, lon = c(-2, -0.5), lat = c(43.2, 43.7))
#' nc2 <- downscaleNcdf(nc, year = 2005, month = 3:8, lon = c(-2, -0.5), lat = c(43.2, 43.7))
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export 
#' @references 
#' 
#' \itemize{
#' 
#' \item Santander MetGroup (2015). ecomsUDG.Raccess: R interface to the ECOMS User Data Gateway. R package
#' version 2.2-6. http://meteo.unican.es/ecoms-udg
#' }
#' 
#' 
downscaleNcdf <- function(gridData, year = NULL, month = NULL, lon = NULL, lat = NULL) {
  
  
  if (!is.null(year)) {
    Dates <- as.POSIXlt(gridData$Dates$start)
    yearIndex <- Dates$year + 1900
    monIndex <- Dates$mon + 1
    timeDim <- match('time', attributes(gridData$Data)$dimensions)
    
    
    if (is.null(month) || !any(sort(month) != month)) {
      targetYearIndex <- which(yearIndex %in% year)
      if (length(targetYearIndex) == 0) stop('No input years in the input ts, check your input.')
      
      
      # if year crossing  than sort(month) != month
    } else {
      
      startIndex <- intersect(which(yearIndex == year[1] - 1), which(monIndex == month[1]))[1]
      endIndex <- tail(intersect(which(yearIndex == tail(year, 1)), which(monIndex == tail(month, 1))), 1)
      
      if (is.na(startIndex) || length(endIndex) == 0 || startIndex > endIndex) {
        stop('Cannot find input months and input years in the input time series.')
      } else {
        
        targetYearIndex <- startIndex:endIndex
        
        if (any(diff(year) != 1)) {
          # if year is not continuous, like 1999, 2003, 2005, than we have to sift again.
          # Only for special cases.
          Dates <- Dates[targetYearIndex]
          yea <- Dates$year + 1900
          mon <- Dates$mon + 1
          
          DateIndex <- unlist(sapply(year, function(x) {
            startIndex <- intersect(which(yea == x - 1), which(mon == month[1]))[1]
            endIndex <- tail(intersect(which(yea == x), which(mon == tail(month, 1))), 1)
            index <- startIndex:endIndex
            return(index)
          }))
          
          
          targetYearIndex <- targetYearIndex[DateIndex]
          # cannot directly return output here, because sometimes, month can be incontinuous,
          # we still need the next process to sift month.
        }
      }
    }
    
    gridData$Dates$start <- gridData$Dates$start[targetYearIndex]
    gridData$Dates$end <- gridData$Dates$end[targetYearIndex]
    
    gridData$Data <- chooseDim(gridData$Data, timeDim, targetYearIndex)
  }  
    
  if (!is.null(month)) {
    Dates <- as.POSIXlt(gridData$Dates$start)
    monIndex <- Dates$mon + 1
    
    targetMonIndex <- which(monIndex %in% month)
    if (length(targetMonIndex) == 0) stop('Check your input year, it may exceed the years 
                                          in the input dataset.')
    gridData$Dates$start <- gridData$Dates$start[targetMonIndex]
    gridData$Dates$end <- gridData$Dates$end[targetMonIndex]
    
    timeDim <- match('time', attributes(gridData$Data)$dimensions)
    
    gridData$Data <- chooseDim(gridData$Data, timeDim, targetMonIndex)
    
  }
  
  if (!is.null(lon)) {
    
    lonIndex <- gridData$xyCoords$x
    
    lonI1 <- which(abs(lonIndex - min(lon)) == min(abs(lonIndex - min(lon)), na.rm = TRUE)) 
    lonI2 <- which(abs(lonIndex - max(lon)) == min(abs(lonIndex - max(lon)), na.rm = TRUE)) 
    
    # take the as large as possible range
    targetLonIndex <- lonI1[length(lonI1)]:lonI2[length(lonI2)]
    if (length(targetLonIndex) == 0) stop('Your input lon is too small, try to expand the 
                                          longitude range.') 
    gridData$xyCoords$x <- gridData$xyCoords$x[targetLonIndex]
    lonDim <- grepAndMatch('lon', attributes(gridData$Data)$dimensions)
    
    gridData$Data <- chooseDim(gridData$Data, lonDim, targetLonIndex)
  }
  
  
  if (!is.null(lat)) {
    latIndex <- gridData$xyCoords$y
    
    latI1 <- which(abs(latIndex - min(lat)) == min(abs(latIndex - min(lat)), na.rm = TRUE)) 
    latI2 <- which(abs(latIndex - max(lat)) == min(abs(latIndex - max(lat)), na.rm = TRUE)) 
    
    targetLatIndex <- latI1[length(latI1)]:latI2[length(latI2)]
    
    if (length(targetLonIndex) == 0) stop('Your input lat is too small, try to expand the 
                                          latitude range.') 
    gridData$xyCoords$y <- gridData$xyCoords$y[targetLatIndex]
    latDim <- grepAndMatch('lat', attributes(gridData$Data)$dimensions)
    gridData$Data <- chooseDim(gridData$Data, latDim, targetLatIndex)
  }
  
  return(gridData)
  
}










#' Write to NetCDF file using hyfo list file
#' @param gridData A hyfo list file from \code{\link{loadNcdf}}
#' @param filePath A path of the new NetCDF file, should end with ".nc"
#' @param missingValue A number representing the missing value in the NetCDF file, default
#' is 1e20
#' #' @param tz A string representing the time zone, default is GMT, if you know what time zone is 
#' you can assign it in the argument. If \code{tz = ''}, current time zone will be taken.
#' @param units A string showing in which unit you are putting in the NetCDF file, it can be 
#' seconds or days and so on. If not specified, the function will pick up the possible largest 
#' time units from \code{c('weeks', 'days', 'hours', 'mins', 'secs')}
#' @param version ncdf file versions, default is 3, if 4 is chosen, output file will be foreced to version 4.
#' @return An NetCDF version 3 file.
#' @examples 
#' # First open the test NETcDF file.
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")
#' 
#' 
#' # Then if you don't know the variable name, you can use \code{getNcdfVar} to get variable name
#' varname <- getNcdfVar(filePath)
#' 
#' nc <- loadNcdf(filePath, varname)
#' 
#' # Then write to your work directory
#' 
#' writeNcdf(nc, 'test.nc')
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export 
#' @import ncdf4
#' @references 
#' 
#' \itemize{
#' \item David Pierce (2015). ncdf4: Interface to Unidata netCDF (Version 4 or
#' Earlier) Format Data Files. R package version 1.14.1.
#' http://CRAN.R-project.org/package=ncdf4
#' 
#' \item Santander MetGroup (2015). ecomsUDG.Raccess: R interface to the ECOMS User Data Gateway. R package
#' version 2.2-6. http://meteo.unican.es/ecoms-udg
#' 
#' }
#' 
#' 
writeNcdf <- function(gridData, filePath, missingValue = 1e20, tz = 'GMT', units = NULL, version = 3) {
  
  name <- gridData$Variable$varName
  # First defines dimensions.
  lonName <- attributes(gridData$xyCoords$x)$name
  latName <- attributes(gridData$xyCoords$y)$name
  dimLon <- ncdim_def(lonName, 'degree', gridData$xyCoords$x)
  dimLat <- ncdim_def(latName, 'degree', gridData$xyCoords$y)
  dimMem <- NULL
  if (!is.null(gridData$Members)) {
    dimMem <- ncdim_def('member', 'members', 1:length(gridData$Members))
  }
  
  
  # Time needs to be treated seperately
  dates <- as.POSIXlt(gridData$Dates$start, tz = tz)
  if (is.null(units)) {
    units <- getTimeUnit(dates)
    time <- difftime(dates, dates[1], units = units)
  } else {
    time <- difftime(dates, dates[1], units = units)
  }
  timeUnits <- paste(units, 'since', dates[1])
  # Here time needs to be numeric, as required by ncdf4 package, which is not the same
  # with ncdf
  dimTime <- ncdim_def('time', timeUnits, as.numeric(time))
  
  
  # Depending on whether there is a member part of the dataset.
  # default list
  dimList <- list(dimLon, dimLat, dimTime, dimMem)
  
  # In order to keep the dim list exactly the same with the original one, it needs to be changed.
  dimIndex <- grepAndMatch(c('lon', 'lat', 'time', 'member'), attributes(gridData$Data)$dimensions)
  dimIndex <- na.omit(dimIndex)
  
  # Here order is needed, cuz in the procesure above, c('lon', 'lat', 'time', 'member')
  # is the pattern, while actually, attributes(gridData$Data)$dimensions should be the pattern.
  # So here needs an order() to get the wanted result.
  dimList <- dimList[order(dimIndex)]
  
  # delete the NULL list, in order that there is no member part in the data.
  dimList <- Filter(Negate(is.null), dimList)
  # Then difines data
  var <- ncvar_def( name, "units", dimList, missingValue)
  
  
  # Here for ncdf4, there is an option to create version 4 ncdf, in future, it
  # may added here.
  if (version == 3) {
    nc <- nc_create(filePath, var) 
  } else if (version == 4) {
    nc <- nc_create(filePath, var, force_v4 = TRUE)
  } else {
    stop("Which ncdf version you want? Only 3 and 4 can be selected!")
  }
  
  # This part comes from the library downscaleR, can be deleted if you don't 
  # use {ecomsUDG.Raccess}, by adding this, the file can be read by the package {ecomsUDG.Raccess}
  ncatt_put(nc, "time", "standard_name","time")
  ncatt_put(nc, "time", "axis","T")
  ncatt_put(nc, "time", "_CoordinateAxisType","Time")
  #ncatt_put(nc, "time", "_ChunkSize",1)
  ncatt_put(nc, lonName, "standard_name","longitude")
  ncatt_put(nc, lonName, "_CoordinateAxisType","Lon")
  ncatt_put(nc, latName, "standard_name","latitude")
  ncatt_put(nc, latName, "_CoordinateAxisType","Lat")
  if (!is.null(dimMem)){
    ncatt_put(nc, "member", "standard_name","realization")
    ncatt_put(nc, "member", "_CoordinateAxisType","Ensemble")
    #att.put.ncdf(nc, "member", "ref","http://www.uncertml.org/samples/realisation")
  }
  
  
  # This part has to be put
  ncatt_put(nc, 0, "Conventions","CF-1.4")
  ncatt_put(nc, 0, 'WrittenBy', 'hyfo(http://yuanchao-xu.github.io/hyfo/)')
  
  #data <- aperm(gridData$Data, dimIndex) no need to do this, in the process above
  # when you define the dimlist, the order of the dimension was fixed.
  data <- gridData$Data
  ncvar_put(nc, name, data)
  nc_close(nc)
  
}

# For internaluse by writeNcdf
getTimeUnit <- function(dates) {
  units <- c('weeks', 'days', 'hours', 'mins', 'secs')
  output <- NULL
  for (unit in units) {
    time <- difftime(dates, dates[1], units = unit)
    rem <- sapply(time, function(x) x%%1)
    if (!any(rem != 0)) {
      output <- unit
      break
    }
  } 
  return(output)
}


# Save for future use. 
#' @import ncdf4
#' @references 
#' David Pierce (2014). ncdf: Interface to Unidata netCDF data files. R package version 1.6.8.
#' http://CRAN.R-project.org/package=ncdf
getExtralDim <- function(...) {
  dimList <- list(...)
  
  
}

# in order to first grep than match.
# match only provides for exactly match, 
# dimIndex <- grepAndMatch(c('lon', 'lat', 'time', 'member'), dimNames)
grepAndMatch <- function(x, table) {
  index <- unlist(lapply(x, function(x) {
    a <- grep(x, table)
  }))
  return(index)
}