station_check <- function(stations){
  if(!(all(stations %in% hail::hydra_data$station))){
    stop("Some of the station names you have provided are not valid")
  }
  stations <- data.frame(station = stations, stringsAsFactors = FALSE)
  results <- merge(stations, hail::hydra_data, all.x = TRUE, all.y = FALSE, by = "station")

  return(list(
    stations = results$station,
    urls = results$url
  ))
}

hydra_base <- function(url, station){

  # Construct URL and read in data
  data <- readLines(url)

  # Munge it
  data <- gsub(x = data[12:length(data)], pattern = " {1,}", perl = TRUE, replacement = "\t")
  data <- utils::read.delim(text = data, header = FALSE, stringsAsFactors = FALSE,
                            col.names = c("date", "daily_total", paste("hour", seq(0,23,1), sep = "_")),
                            na.strings = c("-", "", "NA"))

  data$date <- as.Date(data$date, format = "%d-%B-%Y")

  classes <- unlist(lapply(data, class))
  non_int <- which(classes == "character")
  if(length(non_int)){
    suppressWarnings({
      data[,non_int] <- lapply(data[,non_int], as.integer)
    })
  }
  data$station <- station

  # Done
  return(data)
}

#'@title Read HYDRA data
#'
#'@description \code{hail_hydra} reads hourly rainfall data from
#'the HYDRA network of rain gages around Portland, Oregon. This
#'data is preliminary and raw, meaning that it is likely to contain
#'missing or corrupt values; efforts have been made to ensure such values
#'come out as NAs.
#'
#'@param stations a vector of station names you want the information for;
#'these values can be retrieved from the \code{\link{hydra_data}} dataset.
#'By default, the value is "all", meaning data is retrieved for every station.
#'
#'@return a data.frame containing columns representing the date, daily total
#'rainfall inches, station, and hourly breakdowns.
#'
#'@examples
#'# Simple exmaple
#'walmart_data <- hail_hydra("Walmart Eco Roof")
#'
#'@export
hail_hydra <- function(stations = "all"){

  # If we want all, we want all
  if(stations == "all"){
    stations <- hail::hydra_data$station
  }

  # Check stations are valid
  station_data <- station_check(stations)

  if(length(station_data$urls) > 1){
    results <- mapply(hydra_base, station_data$urls, station_data$stations, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    return(do.call("rbind", results))
  }
  return(hydra_base(station_data$urls, station_data$stations))
}