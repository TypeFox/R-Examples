#' Change lon lat coordinates to cell coordinates
#' @param coord input lon lat coordinate
#' @param lon dataset lon array
#' @param lat dataset lat array
#' @return A cell coordinate
coord2cell <- function(coord, lon, lat) {
  if (length(coord) == 2) {
    # input lon and lat coordinates
    lonC <- coord[1]
    latC <- coord[2]
    
    # Index
    lonI <- which(abs(lon - lonC) == min(abs(lon - lonC), na.rm = TRUE)) 
    latI <- which(abs(lat - latC) == min(abs(lat - latC), na.rm = TRUE))
    
    cell <- c(max(lonI), max(latI))
    
  } else stop('Wrong coord input, should be c(lon, lat). Lon and lat should be within the dataset range.')
  
  return(cell)
}