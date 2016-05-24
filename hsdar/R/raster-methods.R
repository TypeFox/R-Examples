setMethod('extract', signature(x='Speclib', y='ANY'), 
          function(x, y, ...)
{ 
  if (!x@spectra@fromRaster)
    stop("x does not contain RasterBrick")
  
  returnSpeclib <- FALSE
  gotID <- FALSE
  
  res <- extract(x@spectra@spectra_ra, y, ...)
  if (is.list(res))
  {
    if (length(res) > 1)
    {
      n <- 1:length(res)
      n[1] <- nrow(res[[1]])
      for (i in 2:length(res))
        n[i] <- nrow(res[[i]]) + n[i-1]
    } else {
      n <- nrow(res[[1]])
    }
    spec <- matrix(0, ncol = nbands(x), nrow = n[length(n)])

    id <- c(1:nrow(spec))
    
    n <- matrix(c(1:length(n), 1, n[-length(n)] + 1, n ), ncol = 3) 

    res <- apply(n, 1, function(i, x) 
    {
      spec[c(i[2]:i[3]),] <<- as.matrix(x[[i[1]]])
      id[c(i[2]:i[3])] <<- rep.int(i[1], length(i[2]:i[3])) 
    }, res)
    res <- spec
    returnSpeclib <- TRUE
    gotID <- TRUE
  }
  
  if (class(res) == "RasterBrick")
    returnSpeclib <- TRUE
  
  if (class(res) == "data.frame")
  {
    res <- as.matrix(res)
    returnSpeclib <- TRUE
  }
  
  if (returnSpeclib)
  {
    res <- speclib(res, wavelength(x))
    if (gotID)
      attribute(res) <- data.frame(ID = id)
    usagehistory(res) <- usagehistory(x)
    usagehistory(res) <- paste("Values extracted using object of class '",
                               class(y), "' as overlay", sep = "")
  }
  
  return(res)
})
  
setMethod('writeRaster', signature(x='Speclib', filename='character'), 
          function(x, filename, ...) 
{
  if (!x@spectra@fromRaster)
    stop("x does not contain RasterBrick")
  
  spectra(x) <- writeRaster(x@spectra@spectra_ra, filename = filename, ...)
  return(x)
})