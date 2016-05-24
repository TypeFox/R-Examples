#' @title Read DAS File
#' @description Reads a DAS file into a data.frame where each line is data for a specific
#'   event.
#'
#' @param file filename of a DAS file.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
das.read <- function(file) {
  # read data from file as text lines
  DAS <- readLines(file, warn = FALSE)
  nDAS <- length(DAS)

  # read event and effort
  Event <- substr(DAS, 4, 4)
  OnEffort <- (substr(DAS, 5, 5) == ".")

  # date and time
  tm <- substr(DAS, 6, 11)
  tm <- gsub(" ", NA, tm)
  dt <- substr(DAS, 13, 18)
  dt <- gsub(" ", NA, dt)
  Date <- strptime(paste(dt, tm), "%m%d%y %H%M%S")

  # latitude and longitude
  suppressWarnings({
    LatD <- as.numeric(substr(DAS, 21, 22))
    LatM <- as.numeric(substr(DAS, 24, 28))
    LongD <- as.numeric(substr(DAS, 31, 33))
    LongM <- as.numeric(substr(DAS, 35, 39))
  })
  Lat <- (LatD + LatM/60) * ifelse(substr(DAS, 20, 20) == "S", -1, 1)
  Long <- (LongD + LongM/60) * ifelse(substr(DAS, 30, 30) == "W", -1, 1)

  # OTHER EVENT-SPECIFIC DATA FIELDS
  Data1 <- substr(DAS, 40, 44)
  Data2 <- substr(DAS, 45, 49)
  Data3 <- substr(DAS, 50, 54)
  Data4 <- substr(DAS, 55, 59)
  Data5 <- substr(DAS, 60, 64)
  Data6 <- substr(DAS, 65, 69)
  Data7 <- substr(DAS, 70, 74)

  # read species
  Spp1 <- Spp2 <- Spp3 <- rep('', nDAS)
  event.A <- Event == "A"
  Spp1[event.A] <- substr(Data5[event.A], 3, 5)
  Spp1[Spp1 == ""] <- NA
  Spp1 <- gsub("[[:blank:]]", NA, Spp1)
  Spp2[event.A] <- substr(Data6[event.A], 3, 5)
  Spp2[Spp2 == ""] <- NA
  Spp2 <- gsub("[[:blank:]]", NA, Spp2)
  Spp3[event.A] <- substr(Data7[event.A], 3, 5)
  Spp3[Spp3 == ""] <- NA
  Spp3 <- gsub("[[:blank:]]", NA, Spp3)

  # ASSIGN BEAUFORT SEA STATES
  Bft <- rep(NA, nDAS)
  event.V <- Event == "V"
  Bft[event.V] <- as.numeric(substr(Data1[event.V], 3, 5))
  LastBft <- NA
  for(i in 1:nDAS) {
    if(is.na(Bft[i])) Bft[i] <- LastBft else LastBft <- Bft[i]
  }

  # Sighting, Angle, and Distance
  event.S <- Event == "S"
  Sight <- Angle <- Dist <- rep(NA, nDAS)
  Sight[event.S | event.A] <- as.numeric(Data1[event.S | event.A])
  Angle[event.S] <- as.numeric(Data5[event.S])
  Dist[event.S] <- as.numeric(Data7[event.S])

  data.frame(Event, OnEffort, Date, Lat, Long, Spp1, Spp2, Spp3, Bft, Sight, Angle, Dist, stringsAsFactors = FALSE)
}
