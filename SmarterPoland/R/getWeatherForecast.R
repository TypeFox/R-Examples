addWeatherVariables <- function(df) {
  if ("temperature" %in% colnames(df))
    df$temperatureCelsius <- 5*(df$temperature-32)/9
  if ("apparentTemperature" %in% colnames(df))
    df$apparentTemperatureCelsius <- 5*(df$apparentTemperature-32)/9
  if ("time" %in% colnames(df))
    df$time <-  as.POSIXct(df$time, origin = "1970-01-01", tz = df$timezone)  
  df
}
addWeatherDailyVariables <- function(df) {
  if ("temperatureMax" %in% colnames(df))
    df$temperatureMaxCelsius <- 5*(df$temperatureMax-32)/9
  if ("temperatureMin" %in% colnames(df))
    df$temperatureMinCelsius <- 5*(df$temperatureMin-32)/9
  if ("apparentTemperatureMax" %in% colnames(df))
    df$apparentTemperatureMaxCelsius <- 5*(df$apparentTemperatureMax-32)/9
  if ("apparentTemperatureMin" %in% colnames(df))
    df$apparentTemperatureMinCelsius <- 5*(df$apparentTemperatureMin-32)/9
  if ("temperatureMinTime" %in% colnames(df))
    df$temperatureMinTime <-  as.POSIXct(df$temperatureMinTime, origin = "1970-01-01", tz = df$timezone)  
  if ("temperatureMaxTime" %in% colnames(df))
    df$temperatureMaxTime <-  as.POSIXct(df$temperatureMaxTime, origin = "1970-01-01", tz = df$timezone)  
  if ("apparentTemperatureMin" %in% colnames(df))
    df$apparentTemperatureMin <-  as.POSIXct(df$apparentTemperatureMin, origin = "1970-01-01", tz = df$timezone)  
  if ("apparentTemperatureMax" %in% colnames(df))
    df$apparentTemperatureMax <-  as.POSIXct(df$apparentTemperatureMax, origin = "1970-01-01", tz = df$timezone)  
  if ("sunriseTime" %in% colnames(df))
    df$sunriseTime <-  as.POSIXct(df$sunriseTime, origin = "1970-01-01", tz = df$timezone)  
  if ("sunsetTime" %in% colnames(df))
    df$sunsetTime <-  as.POSIXct(df$sunsetTime, origin = "1970-01-01", tz = df$timezone)  
  if ("time" %in% colnames(df))
    df$time <-  as.POSIXct(df$time, origin = "1970-01-01", tz = df$timezone)  
  df
}


getWeatherForecast <- function(apiKey, lat = NA, lon = NA, city = NA, raw = FALSE) {
  if (!is.na(city)) {
    if (!any(city == rownames(SmarterPoland::cities_lon_lat)))
      simpleError("Cannot find coordinates for this city [this option is only for cities >50k]")
    cityInfo <- SmarterPoland::cities_lon_lat[city,]
    forecast <- GET(paste0("https://api.forecast.io/forecast/",apiKey,"/",cityInfo$lat, ",", cityInfo$long))
  } else {
    if (!is.na(lat) & !is.na(lon)) {
      forecast <- GET(paste0("https://api.forecast.io/forecast/",apiKey,"/",lat, ",", lon))
    } else {
      simpleError("You have to specify city or lat/lon")
    }
  }
  
  forecastJson <- rjson::fromJSON(rawToChar(forecast$content), method = "C")
  if (raw) return(forecastJson)
  
  now = addWeatherVariables(as.data.frame(forecastJson$currently))
  if (class(forecastJson$hourly$data) == "list") {
    props <- c("time", "summary", "precipProbability", "temperature", "apparentTemperature", "dewPoint", "humidity", "windSpeed")
    tmp <- do.call(data.frame, lapply(props, function(p) {
      do.call(rbind, lapply(forecastJson$hourly$data, function(x) x[[p]] ))
            }))
    colnames(tmp) = props
    by.hour = addWeatherVariables(tmp)
  } else {
    by.hour = addWeatherVariables(forecastJson$hourly$data)
  }
  if (class(forecastJson$hourly$data) == "list") {
    props <- c("time", "summary", "sunriseTime", "sunsetTime", "temperatureMin", "temperatureMax",
               "apparentTemperatureMax", "apparentTemperatureMin",
               "temperatureMaxTime","temperatureMinTime",
               "dewPoint", "humidity", "windSpeed")
    tmp <- do.call(data.frame, lapply(props, function(p) {
      do.call(rbind, lapply(forecastJson$daily$data, function(x) x[[p]] ))
    }))
    colnames(tmp) = props
    by.day= addWeatherDailyVariables(tmp)
  } else {
    by.day= addWeatherDailyVariables(forecastJson$daily$data)
  }
  
  list(now = now, by.hour = by.hour, by.day = by.day)  
}

