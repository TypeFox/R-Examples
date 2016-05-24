plotArea <- function(area, col) {
  if (!is.list(area)) stop("'area' must be a ist of 4 elements as produced by 'get.area()'")
  
  col <- c(rgb(0, 0, 0, 0), col)
  
  image(area$Lon, area$Lat, area$Area, col = col, add = TRUE)
}