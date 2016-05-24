isBinary <-
function(arg) {
  binary <- TRUE
  arg <- as.matrix(na.exclude(arg))
  if (is.matrix(arg)) {
    u_mat <- unique(arg)
    for (each in u_mat) { 
      if (each != 0 && each != 1) { 
        binary <- FALSE; break
        } 
      } #breaks and returns FALSE as soon as it hits a non-missing, non (0,1) value
    return(binary) 
  } 
  else { 
    print('Blau dimension columns must be coercable to matrix')
  }
}
