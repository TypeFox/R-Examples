flowlength <- function(DEM,outlet = NA) {

  ## data preparation

  stopifnot(is(DEM,"matrix"))

  if(min(as.vector(DEM[!is.na(DEM)])) > 999998) {
    stop("DEM contains unrealistic values (> 999998)")
  }

  if(!is.na(outlet[1]) && !is(outlet,"numeric") && length(outlet) != 2){
    stop("Outlet should be a vector of length 2")     
  }

  DEM[is.na(DEM)] <- 999999
  
  nrow <- dim(DEM)[1]
  ncol <- dim(DEM)[2]

  if(is.na(outlet[1])) {
    outlet[1] <- -1
    outlet[2] <- -1
  } else {
    if((outlet[1] > nrow) || (outlet[1] < 1) || (outlet[2] > ncol) || (outlet[2] < 1)){
      stop("Outlet should represent coordinates c(row, column) in the DEM")
    }
  }

  ## calling the function

  result <- .C("flowlength",
               PACKAGE = "topmodel",
               as.double(DEM),
               result = double(nrow*ncol),
               as.integer(nrow),
               as.integer(ncol),
               as.integer(outlet[1]),
               as.integer(outlet[2]))$result

  ## formatting of the results

  result[result < 0] <- NA

  return(matrix(result, nrow=nrow))

}
