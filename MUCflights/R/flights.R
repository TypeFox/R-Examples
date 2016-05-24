
flights <- function(from = NULL, to = NULL, 
                    path = system.file("MUCflights.RData", 
                                       package = "MUCflights")) {
  
  # Check
  ### Datenbank anlegen
  # m <- dbDriver("SQLite")
  # con <- dbConnect(m, dbname = path)

  # ret <- dbGetQuery(con, "SELECT * FROM MUCflights")
  # dbDisconnect(con)
  load(path)

  # Datumswerte als solche spezifizieren
  ret$stt <- as.POSIXct(strftime(ret$stt, format = "%Y-%m-%d %H:%M:%S.0"))
  ret$ett <- as.POSIXct(strftime(ret$ett, format = "%Y-%m-%d %H:%M:%S.0"))

  from <- as.POSIXct(from)
  to <- as.POSIXct(to)
  ret <- subset(ret, ret$stt >= from & ret$stt <= to)

  # Klasse festlegen
  class(ret) <- c("flights", class(ret))
  return(ret)
}
