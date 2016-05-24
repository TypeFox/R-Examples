bsFormu <- function(name.y, name.x, intercept = TRUE, ...)
{
  if (!is.character(name.y) | !is.character(name.x)) {
    stop("Names must be character vectors.\n")
  } 
  iy <- length(name.y)
  ix <- length(name.x) 
  st <- ifelse(intercept == TRUE, " ~ 1 + ", " ~ 0 + ")

  if (iy == 1) {
    nn <- paste(name.y, st, name.x[1], sep = "")
    if (ix > 1) {
      for (j in 2:ix) {
        nn <- paste(nn, " + ", name.x[j], sep = "")
      }
    }
    formu <- as.formula(nn)
  }
  
  if (iy > 1) {  
    formu <- list()   
    for (i in 1:iy) {
      nn <- paste(name.y[i], st, name.x[1], sep = "")
      if (ix > 1) {
        for (j in 2:ix) {
          nn <- paste(nn, " + ", name.x[j], sep = "")
        }
      }
      formu[[i]] <- as.formula(nn)
    }
  }
  return(formu)
}