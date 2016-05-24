as.SpatialLines <- function(x, ...){
  UseMethod("as.SpatialLines")
}
as.SpatialLines.SpatialStreamNetwork <- function(x, ...) {
  SpatialLines(x@lines, proj4string = x@proj4string)
}

as.SpatialLinesDataFrame.SpatialStreamNetwork <- function(x, ...) {
    SpatialLinesDataFrame(x@lines, x@data, proj4string = x@proj4string)
}



as.SpatialPoints <- function(x, ...){

  UseMethod("as.SpatialPoints")
}

as.SpatialPoints.SpatialStreamNetwork <- function(x, data = "Obs", ...) {
  if(data == "Obs") {
  return(SpatialPoints(x@obspoints@SSNPoints[[1]]@point.coords,
        proj4string = x@proj4string))
  } else {
     ind <- x@predpoints@ID %in% data
  if(sum(ind) ==0) {
       stop(paste0(data, " is not present in STSN"))}
     if(sum(ind) > 1) {
       stop(paste0(data, " resides in more than one slot in the STSN"))}
     j <- which(ind)
     return(SpatialPoints(x@predpoints@SSNPoints[[j]]@point.coords,
        proj4string = x@proj4string))
  }
  return(invisible())
}

as.SpatialPointsDataFrame <- function(x, ...){
  UseMethod("as.SpatialPointsDataFrame")
}

as.SpatialPointsDataFrame.SpatialStreamNetwork <- function(x, data = "Obs", ...) {
  if(data == "Obs") {
        return(SpatialPointsDataFrame(x@obspoints@SSNPoints[[1]]@point.coords,
        x@obspoints@SSNPoints[[1]]@point.data,
        proj4string = x@proj4string))
  } else {
     ind <- x@predpoints@ID %in% data
     if(sum(ind) ==0) {
       stop(paste0(data, " is not present in SSN"))}
     if(sum(ind) > 1) {
       stop(paste0(data, " resides in more than one slot in the SSN"))}
     j <- which(ind)
     return(SpatialPointsDataFrame(x@predpoints@SSNPoints[[j]]@point.coords,
        x@predpoints@SSNPoints[[j]]@point.data,
        proj4string = x@proj4string))
  }
  return(invisible())
}

as.SpatialLinesDataFrame <- function(x, ...) {
    UseMethod("as.SpatialLinesDataFrame")
}

as.SpatialLinesDataFrame.SpatialStreamNetwork <- function(x, ...) {
    sl <- SpatialLines(x@lines, proj4string = x@proj4string)
    data <- x@data
    sldf <- SpatialLinesDataFrame(sl, data, match.ID = FALSE)
    sldf
}
