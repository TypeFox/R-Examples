#' SensusR:  Sensus Analytics
#'
#' Provides access and analytic functions for Sensus data. More information can be found at the
#' following URL:
#' 
#'     https://github.com/MatthewGerber/sensus/wiki
#' 
#' @section SensusR functions:
#' The SensusR functions handle reading, cleaning, plotting, and otherwise analyzing data collected
#' via the Sensus system.
#'
#' @docType package
#' 
#' @name SensusR
NULL

#' Read JSON-formatted Sensus data.
#' 
#' @param path Path to JSON file.
#' @param convert.to.local.timezone Whether or not to convert timestamps to the local timezone.
#' @return All data, listed by type.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
read.sensus.json = function(path, convert.to.local.timezone = TRUE)
{
  local.timezone = Sys.timezone()
  
  # read all lines, only retaining non-empty lines
  con = file(path, open="r")
  lines = as.matrix(readLines(con))
  lines = apply(lines, 1, trim)
  lines = as.matrix(lines[sapply(lines, nchar) > 0])
  close(con)

  # parse each line to json
  lines = apply(lines, 1, function(line)
  {
    json = jsonlite::fromJSON(line)
    
    # set short version of type
    datum.type = strsplit(json$"$type", ",")[[1]][1]
    datum.type = tail(strsplit(datum.type, "[.]")[[1]], n=1)
    json$Type = datum.type
    
    # we no longer need the $type column
    json = json[-which(names(json) %in% c("$type"))]
    
    return(as.data.frame(json, stringsAsFactors = FALSE))
  })
  
  # split up data by type
  types = as.factor(sapply(lines, function(line) { return(line$Type) }))
  data = split(lines, types)
  
  # unlist everything
  for(datum.type in levels(types))
  {
    first.row = data[[datum.type]][[1]]
    column.names = names(first.row)
    
    # build new dataframe for the current data type
    new.data = data.frame(matrix(nrow=length(data[[datum.type]]), ncol=0))
    for(col in column.names)
    {
      col.data = unlist(sapply(data[[datum.type]], function(row,col) { return(row[[col]])}, col))
      new.data[[col]] = col.data
    }
    
    # parse/convert all time stamps
    new.data$Timestamp = strptime(new.data$Timestamp, format = "%Y-%m-%dT%H:%M:%OS", tz="UTC")    
    if(convert.to.local.timezone)
    {
      new.data$Timestamp = lubridate::with_tz(new.data$Timestamp, local.timezone)
    }
    
    # don't need type anymore, since we've group by type
    new.data$Type = NULL
    
    # order by timestamp
    new.data = new.data[order(new.data$Timestamp),]
    
    # filter redundant data by Id and remove Id column
    new.data = new.data[!duplicated(new.data$Id),]
    new.data$Id = NULL
    
    # set class to the datum type, in order for generic plot functions to work
    class(new.data) = c(datum.type, class(new.data))
    
    data[[datum.type]] = new.data
  }
  
  return(data)
}

#' Plot accelerometer data.
#' 
#' @method plot AccelerometerDatum
#' @param x Accelerometer data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$AccelerometerDatum)
plot.AccelerometerDatum = function(x, pch = ".", type = "l", ...)
{ 
  par(mfrow=c(2,2))
  plot.default(x$Timestamp, x$X, main = "Accelerometer", xlab = "Time", ylab = "X", pch = pch, type = type)
  plot.default(x$Timestamp, x$Y, main = "Accelerometer", xlab = "Time", ylab = "Y", pch = pch, type = type)
  plot.default(x$Timestamp, x$Y, main = "Accelerometer", xlab = "Time", ylab = "Z", pch = pch, type = type)
  par(mfrow=c(1,1))
}

#' Plot altitude data.
#' 
#' @method plot AltitudeDatum
#' @param x Altitude data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$AltitudeDatum)
plot.AltitudeDatum = function(x, pch = ".", type = "l", ...)
{
  plot.default(x$Timestamp, x$Altitude, main = "Altitude", xlab = "Time", ylab = "Meters", pch = pch, type = type, ...)
}

#' Plot battery data.
#' 
#' @method plot BatteryDatum
#' @param x Battery data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$BatteryDatum)
plot.BatteryDatum = function(x, pch = ".", type = "l", ...)
{
  plot(x$Timestamp, x$Level, main = "Battery", xlab = "Time", ylab = "Level (%)", pch = pch, type = type, ...)
}

#' Plot cell tower data.
#' 
#' @method plot CellTowerDatum
#' @param x Cell tower data.
#' @param ... Other plotting arguments.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$CellTowerDatum)
plot.CellTowerDatum = function(x, ...)
{
  freqs = plyr::count(x$CellTower)
  if(nrow(freqs) > 0)
  {
    pie(freqs$freq, freqs$x, main = "Cell Tower", ...)
  }
}

#' Plot compass data.
#' 
#' @method plot CompassDatum
#' @param x Compass data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$CompassDatum)
plot.CompassDatum = function(x, pch = ".", type = "l", ...)
{
  plot(x$Timestamp, x$Heading, main = "Compass", xlab = "Time", ylab = "Heading", pch = pch, type = type, ...)
}

#' Plot light data.
#' 
#' @method plot LightDatum
#' @param x Light data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$LightDatum)
plot.LightDatum = function(x, pch = ".", type = "l", ...)
{
  plot(x$Timestamp, x$Brightness, main = "Light", xlab = "Time", ylab = "Level", pch = pch, type = type, ...)
}

#' Plot location data.
#' 
#' @method plot LocationDatum
#' @param x Location data.
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$LocationDatum)
plot.LocationDatum = function(x, ...)
{
  lon = x$Longitude
  lat = x$Latitude
  newmap = rworldmap::getMap(resolution = "high")
  plot(newmap, xlim = range(lon), ylim = range(lat), asp = 1, ...)
  points(lon, lat, col = "red", cex = .6, ...)
}

#' Plot running apps data.
#' 
#' @method plot RunningAppsDatum
#' @param x Apps data.
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$RunningAppsDatum)
plot.RunningAppsDatum = function(x, ...)
{
  freqs = plyr::count(x$Name)
  if(nrow(freqs) > 0)
  {
    pie(freqs$freq, freqs$x, main = "Running Apps", ...)
  }
}

#' Plot screen data.
#' 
#' @method plot ScreenDatum
#' @param x Screen data.
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$ScreenDatum)
plot.ScreenDatum = function(x, ...)
{
  plot(x$Timestamp, x$On, main = "Screen", xlab = "Time", ylab = "On/Off", pch=".", type = "l", ...)
}

#' Plot sound data.
#' 
#' @method plot SoundDatum
#' @param x Sound data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$SoundDatum)
plot.SoundDatum = function(x, pch = ".", type = "l", ...)
{
  plot(x$Timestamp, x$Decibels, main = "Sound", xlab = "Time", ylab = "Decibels", pch = pch, type = type, ...)
}

#' Plot speed data.
#' 
#' @method plot SpeedDatum
#' @param x Speed data.
#' @param pch Plotting character.
#' @param type Line type. 
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$SpeedDatum)
plot.SpeedDatum = function(x, pch = ".", type = "l", ...)
{
  plot(x$Timestamp, x$KPH, main = "Speed", xlab = "Time", ylab = "KPH", pch = pch, type = type, ...)
}

#' Plot telephony data.
#' 
#' @method plot TelephonyDatum
#' @param x Telephony data.
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$TelephonyDatum)
plot.TelephonyDatum = function(x, ...)
{
  freqs = plyr::count(x$PhoneNumber[x$PhoneNumber != ""])
  if(nrow(freqs) > 0)
  {
    pie(freqs$freq, freqs$x, main = "Phone Numbers", ...)
  }
}

#' Plot WLAN data.
#' 
#' @method plot WlanDatum
#' @param x WLAN data.
#' @param ... Other plotting parameters.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(data$WlanDatum)
plot.WlanDatum = function(x, ...)
{
  freqs = plyr::count(x$AccessPointBSSID[x$AccessPointBSSID != ""])
  if(nrow(freqs) > 0)
  {
    pie(freqs$freq, freqs$x, main = "WLAN BSSID", ...)
  }
}

#' Get timestamp lags for a Sensus data frame.
#' 
#' @param data Data to plot lags for (e.g., the result of \code{read.sensus.json}).
#' @return List of lags organized by datum type.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' lags = get.all.timestamp.lags(data)
#' plot(lags[["AccelerometerDatum"]])
get.all.timestamp.lags = function(data)
{
  lags = list()
  for(datum.type in names(data))
  {
    if(nrow(data[[datum.type]]) > 1)
    {
      lags[[datum.type]] = hist(as.numeric(diff(data[[datum.type]]$Timestamp)), main = datum.type, xlab = "Lag (Seconds)")
    }
  }
  
  return(lags)
}

#' Get timestamp lags for a Sensus datum.
#' 
#' @param datum One element of a Sensus data frame (e.g., data$CompassDatum).
#' @return List of lags.
#' @examples
#' data = read.sensus.json(system.file("extdata", "example.data.txt", package="SensusR"))
#' plot(get.timestamp.lags(data$AccelerometerDatum))
get.timestamp.lags = function(datum)
{
  lags = NULL
  if(nrow(datum) > 1)
  {
    lags = hist(as.numeric(diff(datum$Timestamp)), xlab = "Lag (Seconds)", main = "Sensus Data")
  } 
  
  return(lags)
}

#' Trim leading white space from a string.
#' 
#' @param x String to trim.
#' @return Result of trimming.
trim.leading = function (x) sub("^\\s+", "", x)

#' Trim trailing white space from a string.
#' 
#' @param x String to trim.
#' @return Result of trimming.
trim.trailing = function (x) sub("\\s+$", "", x)

#' Trim leading and trailing white space from a string.
#' 
#' @param x String to trim.
#' @return Result of trimming.
trim = function (x) gsub("^\\s+|\\s+$", "", x)




