if ( !isGeneric("calcPathRadDOS") ) {
  setGeneric("calcPathRadDOS", function(x, ...)
    standardGeneric("calcPathRadDOS"))
}
#' Compute path radiance based on the dark object method
#'
#' @description
#' Compute an estimated path radiance for all sensor bands, which can then be 
#' used to roughly correct the radiance values for atmospheric scattering. Path 
#' radiance estimation is based on a dark object method.
#'
#' @param x A Satellite object or the value (scaled count) of a dark object in 
#' \code{bnbr} (e.g. minimum raw count of selected raster \code{bnbr}). If x is 
#' a Satellite object, the value is computed using \code{\link{calcDODN}}.
#' @param bnbr Band number for which DNmin is valid.
#' @param band_wls Band wavelengths to be corrected; \code{data.frame} with min 
#' (max) in first (second) column, see details.
#' @param radm Multiplicative coefficient for radiance transformation (i.e. 
#' slope).
#' @param rada Additive coefficient for radiance transformation (i.e. offset).
#' @param szen Sun zenith angle.
#' @param esun Actual (i.e. non-normalized) TOA solar irradiance, e.g. returned 
#' by \code{\link{calcTOAIrradRadRef}}, \code{\link{calcTOAIrradTable}} or 
#' \code{\link{calcTOAIrradModel}}.
#' @param model Model to be used to correct for 1\% scattering (DOS2, DOS4; must 
#' be the same as used by \code{\link{calcAtmosCorr}}).
#' @param esun_method If x is a Satellite object, name of the method to be used 
#' to compute esun using one of \code{\link{calcTOAIrradRadRef}} ("RadRef"), 
#' \code{\link{calcTOAIrradTable}} ("Table") or \code{\link{calcTOAIrradModel}}
#' ("Model")
#' @param dos_adjust Assumed reflection for dark object adjustment; defaults to 0.01.
#' @param scat_coef Scattering coefficient; defaults to -4.0. 
#' @param use_cpp Logical. If \code{TRUE}, C++ functionality (via \strong{Rcpp}) 
#' is enabled, which leads to a considerable reduction of both computation time
#' and memory usage.
#'  
#' @export calcPathRadDOS
#' 
#' @name calcPathRadDOS
#' 
#' @details 
#' If x is a Satellite object, the minimum raw count value (x) is computed using
#' \code{\link{calcDODN}}. If the TOA solar irradiance is not part of the 
#' Satellite object's metadata, it is computed using 
#' \code{\link{calcTOAIrradRadRef}}, \code{\link{calcTOAIrradTable}} or 
#' \code{\link{calcTOAIrradModel}}.
#'  
#' The dark object subtraction approach is based on an approximation 
#' of the atmospheric path radiance (i.e. upwelling radiation which is 
#' scattered into the sensors field of view, aka haze) using the reflectance of a 
#' dark object (i.e. reflectance ~1\%). 
#' 
#' Chavez (1988) proposed a method which uses the dark object reflectance
#' in one band to predict the corresponding path radiances in all other 
#' \code{band_wls}. This is done using a relative radiance model which depends on 
#' the wavelength and overall atmospheric optical thickness (which is estimated 
#' based on the dark object's DN value). This has the advantage that the path 
#' radiance is actually correlated across different sensor \code{band_wls} and 
#' not computed individually for each band using independent dark objects. He 
#' proposed a relative radiance model which follows a wavelength dependent 
#' scattering that ranges from a power of -4 over -2, -1, -0.7 to -0.5 for very 
#' clear over clear, moderate, hazy to very hazy conditions. The relative 
#' factors are computed individually for each 1/1000 wavelength within each band 
#' range and subsequently averaged over the band as proposed by Goslee (2011).
#' 
#' The atmospheric transmittance towards the sensor (Tv) is approximated by 
#' 1.0 (DOS2, Chavez 1996) or Rayleigh scattering (DOS4, Moran et al. 1992)
#' 
#' The atmospheric transmittance from the sun (Tz) is approximated by the 
#' cosine of the sun zenith angle (DOS2, Chavez 1996) or again using Rayleigh
#' scattering (DOS4, Moran et al. 1992).
#' 
#' The downwelling diffuse irradiance is approximated by 0.0 (DOS2, Chavez 1996)
#' or the hemispherical integral of the path radiance (DOS4, Moran et al. 1992).
#' 
#' Equations for the path radiance are taken from Song et al. (2001).
#' 
#' For some sensors, the band wavelengths are already included. See
#' \code{lutInfo()[grepl("_BANDS", names(lutInfo()$META))]} if your sensor is
#' included. To retrieve a sensor, use \code{lutInfo()$<Sensor ID>_BANDS}.
#' 
#' @references Chavez Jr PS (1988) An improved dark-object subtraction technique 
#' for atmospheric scattering correction of multispectral data. Remote Sensing 
#' of Environment 24/3, doi:10.1016/0034-4257(88)90019-3, available online at
#'  \url{http://www.sciencedirect.com/science/article/pii/0034425788900193}.
#'  
#' Chavez Jr PS (1996) Image-based atmospheric corrections revisited and
#' improved. Photogrammetric Engineering and Remote Sensing 62/9,
#' available online at 
#' \url{http://www.asprs.org/PE-RS-Journals-1996/PE-RS-September-1996.html}.
#'  
#' Goslee SC (2011) Analyzing Remote Sensing Data in R: The landsat 
#' Package. Journal of Statistical Software,43/4, 1-25, available online at 
#' \url{http://www.jstatsoft.org/v43/i04/}.
#' 
#' Moran MS, Jackson RD, Slater PN, Teillet PM (1992) Evlauation of simplified
#' procedures for rretrieval of land surface reflectane factors from satellite
#' sensor output.Remote Sensing of Environment 41/2-3, 169-184, 
#' doi:10.1016/0034-4257(92)90076-V, available online at
#' \url{http://www.sciencedirect.com/science/article/pii/003442579290076V}.
#' 
#' Song C, Woodcock CE, Seto KC, Lenney MP, Macomber SA (2001) Classification 
#' and Change Detection Using Landsat TM Data: When and How to Correct 
#' Atmospheric Effects? Remote Sensing of Environment 75/2, 
#' doi:10.1016/S0034-4257(00)00169-3, available online at
#' \url{http://www.sciencedirect.com/science/article/pii/S0034425700001693}
#'
#' If you refer to Sawyer and Stephen 2014, please note that eq. 5 is wrong.
#' 
#' @seealso This function is used by \code{\link{calcAtmosCorr}} to 
#' compute the path radiance.
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' sat <- calcTOAIrradModel(sat)
#' 
#' bcde <- "B002n"
#' 
#' sat <- calcTOAIrradRadRef(sat, normalize = FALSE)
#' 
#' ## performance of base-R
#' system.time(
#'  val1 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
#'                 bnbr = getSatLNBR(sat, bcde),
#'                 band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
#'                                       LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
#'                 radm = getSatRADM(sat, getSatBCDESolar(sat)),
#'                 rada = getSatRADA(sat, getSatBCDESolar(sat)),
#'                 szen = getSatSZEN(sat, getSatBCDESolar(sat)),
#'                 esun = getSatESUN(sat, getSatBCDESolar(sat)),
#'                 model = "DOS2",
#'                 scat_coef = -4, use_cpp = FALSE)
#' )
#' 
#' ## and Rcpp version
#' system.time(
#'  val2 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
#'                 bnbr = getSatLNBR(sat, bcde),
#'                 band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
#'                                       LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
#'                 radm = getSatRADM(sat, getSatBCDESolar(sat)),
#'                 rada = getSatRADA(sat, getSatBCDESolar(sat)),
#'                 szen = getSatSZEN(sat, getSatBCDESolar(sat)),
#'                 esun = getSatESUN(sat, getSatBCDESolar(sat)),
#'                 model = "DOS2",
#'                 scat_coef = -4, use_cpp = TRUE)
#' )
#' 
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return Satellite object with path radiance for each band in the metadata
#' (W m-2 micrometer-1)
#' 
#' @rdname calcPathRadDOS
#'
setMethod("calcPathRadDOS", 
          signature(x = "Satellite"), 
          function(x, model = c("DOS2", "DOS4"), esun_method = "RadRef", 
                   use_cpp = TRUE){

            # if not supplied, model defaults to DOS2
            model <- model[1]
            
            # Compute TOA solar irradiance information if necessary
            if(any(is.na(getSatESUN(x, getSatBCDESolar(x))))){
              # Compute toa irradiance for all solar band layers
              if(esun_method == "Table"){
                x <- calcTOAIrradTable(x, normalize = TRUE)
              } else if(esun_method == "Model"){
                x <- calcTOAIrradModel(x, model = "MNewKur", normalize = TRUE)
              } else if(esun_method == "RadRef"){
                x <- calcTOAIrradRadRef(x, normalize = "TRUE")
              }
            }
            
            # Get solar bands with calibration information equals scaled counts
            sc_bands <- getSatBCDESolarCalib(x, calib = "SC")
            
            # Take care of dark object values
            bcde <- "B002n"
            dn_min <- calcDODN(getSatDataLayer(x, bcde))
            
            # Take care of path radiance
            path_rad <- calcPathRadDOS(x = dn_min,
                                       bnbr = getSatLNBR(x, bcde),
                                       band_wls = data.frame(
                                         LMIN = getSatLMIN(x, sc_bands),
                                         LMAX = getSatLMAX(x, sc_bands)),
                                       radm = getSatRADM(x, sc_bands),
                                       rada = getSatRADA(x, sc_bands),
                                       szen = getSatSZEN(x, sc_bands),
                                       esun = getSatESUN(x, sc_bands),
                                       model = model, 
                                       use_cpp = use_cpp)
            x <- addSatMetaParam(x, 
                                 meta_param = data.frame(
                                   BCDE = names(path_rad),
                                   PRAD = as.numeric(path_rad)))
            return(x)
          })


# Function using numeric -------------------------------------------------------
#' 
#' @return Vector object with path radiance values for each band 
#' (W m-2 micrometer-1)
#' 
#' @rdname calcPathRadDOS
#'
setMethod("calcPathRadDOS", 
          signature(x = "numeric"), 
          function(x, bnbr, band_wls, radm, rada, szen, esun,
                   model = c("DOS2", "DOS4"), 
                   scat_coef = c(-4.0, -2.0, -1.0, -0.7, -0.5), 
                   dos_adjust = 0.01, use_cpp = TRUE){
            
            # if not supplied, model defaults to DOS2
            model <- model[1]
            
            # if not supplied, scat_coef defaults to -4.0
            scat_coef <- scat_coef[1]
            
            # Define relative scattering model based on wavelength dependent 
            # scattering processes and different atmospheric optical thiknesses. 
            # Resulting  coefficients are the band wavelength to the power of 
            # the respective scattering coefficients defined for the model.
            if (use_cpp) {
              band_wls_mat <- as.matrix(band_wls)
              scattering_model <- ScatteringModel(band_wls_mat, scat_coef)
            } else {
              scattering_model <- sapply(seq(nrow(band_wls)), function(y){
                act_band <- seq(band_wls[y, 1], band_wls[y, 2], by=0.001)
                mean(act_band ^ scat_coef)
              })
            }
            
            # Compute multiplication factors which relate the scattering model 
            # coefficents between the individual band_wls and the one band which 
            # has been used for the eytraction of the dark object radiation.
            basline_factor <- scattering_model[bnbr]
            scattering_factors <- scattering_model / basline_factor
            
            # Compute radiation multiplication factors for radiance conversion 
            # (i.e.RADM) and normalize the factors to the multiplication factor 
            # of the band that has been used for the dark object eytraction.
            normalized_radm <- radm / radm[bnbr]
            
            # Compute first estimate of path radiance for all band_wls based on 
            # the values the band used to derive the dark object properties
            lp_min <- radm[bnbr] * x + rada[bnbr]
            lp_min_band_wls <- lp_min * scattering_factors
            
            # Compute a correction term for the path radiance values to consider
            # a minimum surface reflection wich one can always eypect.
            esun
            cos_szen <- cos(szen * pi / 180.0)
            if(model == "DOS2"){
              tv <- 1.0
              tz <- cos_szen
              edown <- 0.0
            } else if(model == "DOS4"){
              tv <- NA
              tz <- NA
              edown <- NA
            }
            lp_adj <- dos_adjust * (esun * cos_szen * tz + edown) * tv / pi
            
            # Compute final path radiance estimate for all band_wls
            lp <- lp_min_band_wls - lp_adj
            return(lp)
          })
          