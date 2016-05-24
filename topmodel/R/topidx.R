topidx <- function(DEM, resolution, river=NULL) {

  ## data preparation

  nrow <- dim(DEM)[1]
  ncol <- dim(DEM)[2]

  ## data checking

  stopifnot(is(DEM, "matrix"),(is(river,"matrix") || is.null(river)))

  if(min(as.vector(DEM[!is.na(DEM)])) < -9000)
     stop("DEM contains unrealistic values (< -9000)")
  
  if(is.null(river)) { river = rep(0,nrow*ncol)
  } else { if(min(river) < 0) stop("Error: the object 'river' should only contain positive values") }

  DEM[is.na(DEM)] <- -9999

  ## calling the function

  result <- .C("topidx",
               PACKAGE = "topmodel",
               as.double(DEM),
               as.integer(river),
               as.integer(nrow),
               as.integer(ncol),
               as.double(resolution),
               as.double(resolution),
               result = double(length(DEM)*2))$result

  ## formatting of the results

  atb  <- matrix(result[1:(nrow*ncol)],nrow=nrow)
  area <- matrix(result[(nrow*ncol+1):(nrow*ncol*2)],nrow=nrow)

  atb[atb < -9] <- NA
  area[area < -9] <- NA

  return(list(atb = atb,area = area))

}
