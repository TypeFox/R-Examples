##' Impute missing data points
##'
##' Replace \code{NA}s in the spectra matrix by linear interpolation.
##' @rdname spc.NA.linapprox
##' @param spc hyperSpec object with spectra matrix containing \code{NA}s
##' @param neighbours how many neighbour data points should be used to fit the line
##' @param ... ignored
##' @return hyperSpec object
##' @export
##' @author Claudia Beleites
spc.NA.linapprox <- function (spc, neighbours = 1, ...){
  chk.hy (spc)
  validObject (spc)
  
  all.na <- which (apply (is.na (spc@data$spc), 1, all))
  if (length (all.na) > 0){
    warning ("Spectra containing only NAs found. They will not be approximated.")
  }
  
  ispc <- which (is.na (spc@data$spc), arr.ind = TRUE)

  ispc <- setdiff (unique (ispc[,"row"]), all.na)

  for (i in ispc){
    nas <- which (is.na (spc@data$spc[i,]))
    
    start <- c (0, which (diff (nas) > 1)) + 1
    end   <- c (start [-1] - 1, length (nas)) 

    for (j in seq (along = start)) {
      pts <- nas [start [j]] : nas [end [j]]

      xneighbours <- c ( -(1 : neighbours) + nas [start [j]],
                          (1 : neighbours) + nas [end   [j]]) 
      xneighbours <- xneighbours [xneighbours > 0]
      xneighbours <- xneighbours [xneighbours < nwl (spc)]

      if (length (xneighbours) == 0) # should not happen as all NA-only spectra were excluded
        stop ("No data to interpolate from.") 
      else if (length (xneighbours) == 1)
        spc@data$spc [i, pts] <- spc@data$spc [i, xneighbours]
      else
        spc@data$spc [i, pts] <- approx (x = spc@wavelength  [xneighbours],
                                         y = spc@data$spc [i, xneighbours],
                                         xout = spc@wavelength  [pts],
                                         method = "linear",
                                         rule = 2)$y
      
    }
  }

  spc
}
