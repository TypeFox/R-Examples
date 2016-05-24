GeoCount <- function (x, y, distance, x.latitude = "latitude",
                      x.longitude = "longitude", y.latitude = "latitude",
                      y.longitude = "longitude") {

  if (!(x.latitude %in% colnames(x)))
    stop(paste("Variable", x.latitude, "does not exist in x"))
  if (!(x.longitude %in% colnames(x)))
    stop(paste("Variable", x.longitude, "does not exist in x"))
  if (!(y.latitude %in% colnames(y)))
    stop(paste("Variable", y.latitude, "does not exist in y"))
  if (!(y.longitude %in% colnames(y)))
    stop(paste("Variable", y.longitude, "does not exist in y"))
  

  n.x <- nrow(x)
  n.y <- nrow(y)

  x.lon <- eval(parse(text = paste("x$", x.longitude, sep ="")))
  x.lat <- eval(parse(text = paste("x$", x.latitude, sep = "")))

  y.lon <- eval(parse(text = paste("y$", y.longitude, sep = "")))
  y.lat <- eval(parse(text = paste("y$", y.latitude, sep = "")))

  
  temp <- .C("R2GeoCount",
             as.numeric(x.lon),
             as.numeric(x.lat),
             as.integer(n.x),
             as.numeric(y.lon),
             as.numeric(y.lat),
             as.integer(n.y),
             as.numeric(distance),
             Store.count = integer(n.x),
             package = "endorse")

  res <- as.double(temp$Store.count)
  
  return(res)
}
