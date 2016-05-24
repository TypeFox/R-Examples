rastermeta <- function(x, dim, ext, crs)
{
  if (any(c(missing(dim),
            missing(ext),
            missing(crs))))
  {
    if (missing(x))
      stop("'x' required if not all other arguments are passed")
    
    if (!(class(x) %in% c("Raster", 
                          "RasterBrick", 
                          "Raster", 
                          "HyperSpecRaster",
                          "RasterStack")))
      stop("'x' does not contain all relevant information")
      
    if (missing(dim))
    {
      dim <- c(nrow(x), ncol(x))
    } else {
      if (length(dim) != 2)
        stop("Invalid number of dimensions")
      if (!is.numeric(dim))
        stop("'dim' must be numeric")
    }
    if (missing(ext))
    {
      ext <- extent(x)
    } else {
      if (class(ext) != "extent")
        stop("'ext' must be object of class 'extent'")
    }
    if (missing(crs))
    {
      crs <- crs(x)
    } else {
      if (class(crs) != "CRS")
        stop("'crs' must be object of class 'CRS'")
    }
  } else {
    if (length(dim) != 2)
      stop("Invalid number of dimensions")
    if (!is.numeric(dim))
      stop("'dim' must be numeric")
    if (class(ext) != "extent")
      stop("'ext' must be object of class 'extent'")
    if (class(crs) != "CRS")
      stop("'crs' must be object of class 'CRS'")    
  }
  return(list(dim = dim,
              ext = ext,
              crs = crs))
}

.is.rastermeta <- function(x)
{
  if (class(x) == "list")
  {
    if (length(x) < 3)
      return(FALSE)
    if (names(x)[1] != "dim")
      return(FALSE)
    if (names(x)[2] != "ext")
      return(FALSE)
    if (names(x)[3] != "crs")
      return(FALSE)
  } else {
    if (class(x) %in% c("Speclib", "Specfeat", "Clman"))
      return(.is.rastermeta(x@rastermeta))
    return(FALSE)
  }
  return(TRUE)
}
      
      