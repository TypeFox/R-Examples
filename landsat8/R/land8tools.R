#' Conversion to TOA Radiance
#' @description Conversion to TOA radiance of satellite data.
#' @param x Image to be converted, in matrix, data frame, or SpatialGridDataFrame format.
#' @param Ml band specific multiplicative rescaling factor from the metadata (MTL file) (RADIANCE_MULT_BAND_x, where x is the band number).
#' @param Al Mp band specific additive rescaling factor from the metadata (MTL file) (RADIANCE_ADD_BAND_x, where x is the band number).
#' @return TOA spectral radiance (Watts/(m2*srad*micro-m)).
#' @author Alexandre dos Santos
#' @references U.S. Geological Survey. 2015.
#'             Landsat 8 (L8) data users handbook. Version 1.0. 97p.
#' @examples
#' data(band5)
#' band5.dn<- as(band5, 'SpatialGridDataFrame')
#' band5.rad<-radconv(band5.dn,5.9150E-03,-29.57525)

radconv<-function(x, Ml, Al)
{
  results <- x
  x <- as.vector(as.matrix(x))
  x<-Ml*x+Al
  if (class(results) == "SpatialGridDataFrame")
    results@data[, 1] <- x
  else if (is.data.frame(x))
    results <- data.frame(matrix(x, nrow = nrow(results),
                                 ncol = ncol(results)))
  else results <- x
  results
}

#' Conversion to TOA Reflectance
#' @description Conversion to TOA reflectance of satellite data.
#' @param x Image to be converted, in matrix, data frame, or SpatialGridDataFrame format.
#' @param Mp Band specific multiplicative rescaling factor from the metadata (MTL file) (REFLECTANCE_MULT_BAND_x, where x is the band number).
#' @param Ap Band specific additive rescaling factor from the metadata (MTL file) (REFLECTANCE_ADD_BAND_x, where x is the band number).
#' @return TOA spectral radiance.
#' @author Alexandre dos Santos
#' @references U.S. Geological Survey. 2015.
#'             Landsat 8 (L8) data users handbook. Version 1.0. 97p.
#' @examples
#' data(band5)
#' band5.dn<- as(band5, 'SpatialGridDataFrame')
#' band5.refl<-reflconv(band5.dn,2.0000E-05,-0.100000)


reflconv<-function(x, Mp, Ap)
{
  results <- x
  x <- as.vector(as.matrix(x))
  x<-Mp*x+Ap
  if (class(results) == "SpatialGridDataFrame")
    results@data[, 1] <- x
  else if (is.data.frame(x))
    results <- data.frame(matrix(x, nrow = nrow(results),
                                 ncol = ncol(results)))
  else results <- x
  results
}

#' Conversion to TOA Reflectance with a Correction for the Sun Angle
#' @description Conversion to TOA reflectance with a correction for the sun angle of satellite data.
#' @param x Image to be converted, in matrix, data frame, or SpatialGridDataFrame format.
#' @param Mp band specific multiplicative rescaling factor from the metadata (MTL file) (REFLECTANCE_MULT_BAND_x, where x is the band number).
#' @param Ap band specific additive rescaling factor from the metadata (MTL file) (REFLECTANCE_ADD_BAND_x, where x is the band number).
#' @param sunelev Sun elevation in degrees is provided in the metadata (MTL file) (SUN_ELEVATION).
#' @return TOA spectral radiance with a correction for the sun angle.
#' @author Alexandre dos Santos
#' @references  U.S. Geological Survey. 2015.
#'             Landsat 8 (L8) data users handbook. Version 1.0. 97p.
#' @examples
#' data(band5)
#' band5.dn<- as(band5, 'SpatialGridDataFrame')
#' band5.reflS<-reflconvS(band5.dn,2.0000E-05,-0.100000,41.12846745)


reflconvS<-function(x, Mp, Ap, sunelev)
{
  results <- x
  x <- as.vector(as.matrix(x))
  suntheta <- cos(90-sunelev)
  x<-(Mp*x+Ap)/suntheta
  if (class(results) == "SpatialGridDataFrame")
    results@data[, 1] <- x
  else if (is.data.frame(x))
    results <- data.frame(matrix(x, nrow = nrow(results),
                                 ncol = ncol(results)))
  else results <- x
  results
}

#' Conversion to At Satellite Brightness Temperature
#' @description Conversion to At satellite brightness temperature of satellite data.
#' @param x Image to be converted, in matrix, data frame, or SpatialGridDataFrame format.
#' @param Ml band specific multiplicative rescaling factor from the metadata (MTL file) (RADIANCE_MULT_BAND_x, where x is the band number).
#' @param Al Mp band specific additive rescaling factor from the metadata (MTL file) (RADIANCE_ADD_BAND_x, where x is the band number).
#' @param K1 band specific thermal conversion constant from the metadata (MTL file) (K1_CONSTANT_BAND_x, where x is the band number, 10 or 11).
#' @param K2 band specific thermal conversion constant from the metadata (MTL file) (K2_CONSTANT_BAND_x, where x is the band number, 10 or 11).

#' @return At satellite brightness temperature in Kelvin (K).
#' @author Alexandre dos Santos
#' @references U.S. Geological Survey. 2015.
#'             Landsat 8 (L8) data users handbook. Version 1.0. 97p.
#' @examples
#' data(band11)
#' band11.dn<- as(band11, 'SpatialGridDataFrame')
#' band11.tempK<-tempconv(band11.dn,3.3420E-04,0.10000, 480.89, 1201.14)


tempconv<-function(x, Ml, Al, K1, K2)
{
  results <- x
  x <- as.vector(as.matrix(x))
  x <-K2 / log(K1/(Ml*x+Al)+1)
  if (class(results) == "SpatialGridDataFrame")
    results@data[, 1] <- x
  else if (is.data.frame(x))
    results <- data.frame(matrix(x, nrow = nrow(results),
                                 ncol = ncol(results)))
  else results <- x
  results
}





