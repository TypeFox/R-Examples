river <- function(DEM,atb,area,res,thatb,tharea) {

  ## data preparation
  stopifnot(is(DEM,"matrix"),is(atb,"matrix"),is(area,"matrix"))

  if(min(as.vector(DEM[!is.na(DEM)])) < -9000)
    stop("DEM contains unrealistic values (< -9000)")

  if(min(as.vector(atb[!is.na(atb)])) < -9000)
    stop("Topographic index map contains unrealistic values (< -9000)")

  if(min(as.vector(area[!is.na(area)])) < -9000)
    stop("Cumulative area map contains unrealistic values (< -9000)")

  DEM[is.na(DEM)] <- -9999
  atb[is.na(atb)] <- -9999
  area[is.na(area)] <- -9999
  
  nrow <- dim(DEM)[1]
  ncol <- dim(DEM)[2]

  ## calling the function
  result <- .C("findrivers",
               PACKAGE = "topmodel",
               as.double(DEM),
               as.double(atb),
               as.double(area),
               result = double(nrow * ncol),
               as.integer(nrow),
               as.integer(ncol),
               as.double(res),
               as.double(thatb),
               as.double(tharea))$result

  ## formatting the results

  result[result < 0] <- NA
  result <- matrix(result, nrow=nrow)

  return(result)

}
