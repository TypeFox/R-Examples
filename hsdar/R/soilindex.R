soilindex <- function(
                      x,
                      index,
                      returnHCR = "auto",
                      weighted = TRUE,
                      ...
             )
{
soilindex_available <- function()
{
  av <- c("SWIR SI", "NSMI", "RI_TM", "RI",
          "NDI", "BI_TM", "SI_TM",
          "HI_TM", "CI_TM"
          )
  return(sort(av))
}  

return_index <- function(x)
{
  if (eval.parent(convertSpatialGrid))
  {
    spec <- speclib(x, 1)
    spec@rastermeta <- gridMeta
    result <- HyperSpecRaster(spec)
  }
  return (x)
}

if (length(names(match.call()))==0)
{
  return(soilindex_available())
}  

if (x@spectra@fromRaster)
  return(.blockwise(speclib_obj =  "x", pos = 1))

x_back <- x

if (!is.speclib(x))
  stop("x is not of class 'Speclib'")
if (!x@continuousdata)
  stop("x does not contain continuous spectra")
if (returnHCR == "auto")
  returnHCR <- .is.rastermeta(x)

convertSpatialGrid <- returnHCR
gridMeta <- x@rastermeta

if (returnHCR)
{
  if (!.is.rastermeta(x))
    stop("If returnHCR, x must contain meta information")
}

if (length(index)>1)
{
  result <- as.data.frame(matrix(data = NA,
                                 nrow = dim(x)[1],
                                 ncol = length(index)))
  for (i in 1:length(index))
  {
    temp <- soilindex(x, index[i], returnHCR=FALSE)
    if (!is.null(temp))
    {
      result[,i] <- temp
    }
  }
  if (nspectra(x) > 1)
  {
    names(result) <- index
    row.names(result) <- idSpeclib(x)
  }
  if (returnHCR)
  {
    spec <- speclib(result, c(1:ncol(result)))
    if (.is.rastermeta(x))
      spec@rastermeta <- x@rastermeta
    result <- HyperSpecRaster(spec)
  }
  return(result)
}

d_indexs <- c()
m <- c(rep.int(1,length(d_indexs)))

# index_current <<- index
# row_names_x <<- row.names(x$spectra)

if (any(index==d_indexs)) 
  x <- derivative.speclib(x, m=m[d_indexs==index], ...)

y <- spectra(x)
x <- wavelength(x)

## Pre-defined indices
if (index=="SWIR SI")
{
  return(return_index(-41.59* (get_reflectance(y,x,2210,weighted)-get_reflectance(y,x,2090,weighted)) +
                1.24*(get_reflectance(y,x,2280,weighted)-get_reflectance(y,x,2090,weighted)) + 0.64))
}

if (index=="RI_TM")
{
  x_TM <- try(spectra(spectralResampling(x_back, "Landsat5")), silent = TRUE)
  if (inherits(x_TM, "try-error"))
  {
    warning("Unable to resample to Landsat 5 for RI_TM calculation")
    return(NULL)
  }
  return(return_index(x_TM[,3]^2/(x_TM[,1]*x_TM[,2]^3)))
}

if (index=="BI_TM")
{
  x_TM <- try(spectra(spectralResampling(x_back, "Landsat5")), silent = TRUE)
  if (inherits(x_TM, "try-error"))
  {
    warning("Unable to resample to Landsat 5 for RI_TM calculation")
    return(NULL)
  }
  return(return_index(((x_TM[,1]^2+x_TM[,2]^2+x_TM[,3]^2)/3)^0.5))
}
if (index=="SI_TM")
{
  x_TM <- try(spectra(spectralResampling(x_back, "Landsat5")), silent = TRUE)
  if (inherits(x_TM, "try-error"))
  {
    warning("Unable to resample to Landsat 5 for RI_TM calculation")
    return(NULL)
  }
  return(return_index((x_TM[,3]-x_TM[,1])/(x_TM[,3]+x_TM[,1])))
}
if (index=="HI_TM")
{
  x_TM <- try(spectra(spectralResampling(x_back, "Landsat5")), silent = TRUE)
  if (inherits(x_TM, "try-error"))
  {
    warning("Unable to resample to Landsat 5 for RI_TM calculation")
    return(NULL)
  }
  return(return_index((2*x_TM[,3]-x_TM[,2]-x_TM[,1])/(x_TM[,2]-x_TM[,1])))
}
if (index=="CI_TM")
{
  x_TM <- try(spectra(spectralResampling(x_back, "Landsat5")), silent = TRUE)
  if (inherits(x_TM, "try-error"))
  {
    warning("Unable to resample to Landsat 5 for RI_TM calculation")
    return(NULL)
  }
  return(return_index((x_TM[,3]-x_TM[,2])/(x_TM[,3]+x_TM[,2])))
}
if (index=="NSMI")
{
  return(return_index((get_reflectance(y,x,1800,weighted) - get_reflectance(y,x,2119,weighted))/
                      (get_reflectance(y,x,1800,weighted) + get_reflectance(y,x,2119,weighted))))
}
if (index=="RI")
{
  return(return_index(get_reflectance(y,x,693,weighted)^2/
                      (get_reflectance(y,x,447,weighted) * get_reflectance(y,x,556,weighted)^3)))
}
if (index=="NDI")
{
  return(return_index((get_reflectance(y,x,840,weighted)-get_reflectance(y,x,1650,weighted))/
                      (get_reflectance(y,x,840,weighted)+get_reflectance(y,x,1650,weighted))))
}

## Self-defining indices
index <- gsub("R", "", gsub("(R[0-9]+)", "get_reflectance(y,x,\\1,weighted)", index, 
                            perl = TRUE)
              )
index <- gsub("D", "", gsub("(D[0-9]+)", "get_reflectance(spectra(derivative.speclib(x_back, m=1, ...)),x,\\1,weighted)", index, 
                            perl = TRUE)
              )
index_val <- try(return_index(eval(parse(text = index))), silent = TRUE)
if (inherits(index_val, "try-error"))
{
  cat("Error in self-defined index string or unimplemented index selected\n")
  cat("Index string evals to:\n")
  cat(paste(index, "\n"))
  return(NULL)
}  
return(index_val)
}