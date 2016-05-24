downloadCANSIM <- function(cansimTableNumber){
  temp <- tempfile()
  url <- "http://www20.statcan.gc.ca/tables-tableaux/cansim/csv/"
  filename <- paste0("0", cansimTableNumber, "-eng")
  url <- paste0(url, filename, ".zip")
  download.file(url, temp, quiet = TRUE)
  data <- read.csv(unz(temp, paste0(filename, ".csv") ), stringsAsFactors = FALSE)
  unlink(temp)

  data <- createStatCanVariables(data)

  data$Vector <- NULL
  data$Coordinate <- NULL
  suppressWarnings(data$Value <- as.numeric(data$Value))

  return(data)
}
