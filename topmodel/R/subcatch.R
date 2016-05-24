subcatch <- function(DEM,outlet) {

  ## data preparation

  stopifnot(is(DEM,"matrix"))

  if(min(as.vector(DEM[!is.na(DEM)])) > 999998)
    stop("DEM contains unrealistic values (> 999999)")

  DEM[is.na(DEM)] <- 999999
  
  nrow <- dim(DEM)[1]
  ncol <- dim(DEM)[2]

  if((outlet[1] > nrow) || (outlet[1] < 1) || (outlet[2] > ncol) || (outlet[2] < 1)){
    stop("Outlet should represent coordinates c(row, column) in the DEM")
  }

  ## calling the function

  result <- .C("subcatch",
               PACKAGE = "topmodel",
               as.double(DEM),
               result = integer(nrow*ncol),
               as.integer(nrow),
               as.integer(ncol),
               as.integer(outlet[1]),
               as.integer(outlet[2]))$result

  ## formatting of the results

  return(matrix(result, nrow=nrow))

}
