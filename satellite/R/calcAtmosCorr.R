if ( !isGeneric("calcAtmosCorr") ) {
  setGeneric("calcAtmosCorr", function(x, ...)
    standardGeneric("calcAtmosCorr"))
}
#' Atmospheric correction of remote sensing data
#'
#' @description
#' The function computes an atmospheric scattering correction and converts
#' the sensors digital numbers to reflectances using
#' \itemize{
#'   \item absolute radiance correction
#'   \item DOS2: a dark object substraction model by Chavez (1996)
#'   \item DOS4: a dark object substratcion model by Moran et al. (1992)
#' }
#'
#' @param x Satellite or Raster* object providing the radiance at the sensor.
#' @param model Model to be used to correct for 1\% scattering (DOS2, DOS4).
#' @param esun_method If x is a Satellite object, name of the method to be used 
#' to compute \code{esun} using one of \code{\link{calcTOAIrradRadRef}} 
#' ("RadRef"), \code{\link{calcTOAIrradTable}} ("Table") or 
#' \code{\link{calcTOAIrradModel}} ("Model").
#' @param path_rad Path radiance, e.g. returned from 
#' \code{\link{calcPathRadDOS}}.
#' @param esun Actual (i.e. non-normalized) TOA solar irradiance, e.g. returned 
#' from \code{\link{calcTOAIrradRadRef}}, \code{\link{calcTOAIrradTable}} 
#' or \code{\link{calcTOAIrradModel}}.
#' @param szen Sun zenith angle.
#'
#' @export calcAtmosCorr
#' 
#' @name calcAtmosCorr
#' 
#' @details 
#' If a Satellite object is passed to the function, and if the required
#' pre-processing has not been performed already, the path radiance is computed 
#' based on a dark object's scaled count value using 
#' \code{\link{calcPathRadDOS}} which will also take care of the TOA solar 
#' irradiance by calling \code{\link{calcTOAIrradModel}}, 
#' \code{\link{calcTOAIrradRadRef}} or \code{\link{calcTOAIrradTable}} 
#' (depending on \code{esun_method}) if necessary. The bands' scaled counts are 
#' converted to radiance using \code{\link{convSC2Rad}}.
#'  
#' The radiometric correction is based on a dark object approach using
#' either the DOS2 (Chavez 1996) or DOS4 (Moran et al. 1992) model.
#' 
#' The minimum reflectance values for the dark object are identified using the
#' approximation of Chavez (1988, see \code{\link{calcPathRadDOS}} for details).
#' 
#' The estimated values of the solar irradiance required for the path radiance
#' can be computed by one of \code{\link{calcTOAIrradTable}} which is used to
#' get readily published values of ESun, \code{\link{calcTOAIrradRadRef}} which 
#' computes ESun based on the actual radiance and reflectance in the scene, or  
#' \code{\link{calcTOAIrradModel}} which computes ESun based on  look-up tables 
#' for the sensor's relative spectral response and solar irradiation spectral data.
#' 
#' The atmospheric transmittance towards the sensor (Tv) is approximated by 
#' 1.0 (DOS2, Chavez 1996) or Rayleigh scattering (DOS4, Moran et al. 1992).
#' 
#' The atmospheric transmittance from the sun (Tz) is approximated by the 
#' cosine of the sun zenith angle (DOS2, Chavez 1996) or again using Rayleigh
#' scattering (DOS4, Moran et al. 1992).
#' 
#' The downwelling diffuse irradiance is approximated by 0.0 (DOS2, Chavez 1996)
#' or the hemispherical integral of the path radiance (DOS4, Moran et al. 1992).
#' 
#' Equations are taken from Song et al. (2001).
#'    
#' @references Chavez Jr PS (1988) An improved dark-object subtraction technique 
#' for atmospheric scattering correction of multispectral data. Remote Sensing 
#' of Environment 24/3, doi:10.1016/0034-4257(88)90019-3, available online at
#'  \url{http://www.sciencedirect.com/science/article/pii/0034425788900193}
#'  
#' Chavez Jr PS (1996) Image-based atmospheric corrections revisited and
#' improved. Photogrammetric Engineering and Remote Sensing 62/9,
#' available online at 
#' \url{http://www.asprs.org/PE-RS-Journals-1996/PE-RS-September-1996.html}
#'  
#' Goslee SC (2011) Analyzing Remote Sensing Data in R: The landsat 
#' Package. Journal of Statistical Software,43/4, 1-25. URL 
#' \url{http://www.jstatsoft.org/v43/i04/}.
#' 
#' Moran MS, Jackson RD, Slater PN, Teillet PM (1992) Evaluation of simplified
#' procedures for retrieval of land surface reflectance factors from satellite
#' sensor output.Remote Sensing of Environment 41/2-3, 169-184, 
#' doi:10.1016/0034-4257(92)90076-V, 
#' URL \url{http://www.sciencedirect.com/science/article/pii/003442579290076V}.
#' 
#' Song C, Woodcock CE, Seto KC, Lenney MP, Macomber SA (2001) Classification 
#' and Change Detection Using Landsat TM Data: When and How to Correct 
#' Atmospheric Effects? Remote Sensing of Environment 75/2, 
#' doi:10.1016/S0034-4257(00)00169-3, URL
#' \url{http://www.sciencedirect.com/science/article/pii/S0034425700001693}
#' 
#'
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' sat_atmos <- calcAtmosCorr(sat, model = "DOS2", esun_method = "RadRef")
#' 
#' bcde <- "B002n"
#' 
#' sat <- calcTOAIrradRadRef(sat, normalize = FALSE)
#' 
#' path_rad <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
#'                            bnbr = getSatLNBR(sat, bcde),
#'                            band_wls = 
#'                              data.frame(LMIN = 
#'                                           getSatLMIN(sat, 
#'                                                      getSatBCDESolar(sat)), 
#'                                         LMAX = 
#'                                           getSatLMAX(sat, 
#'                                                      getSatBCDESolar(sat))),
#'                            radm = getSatRADM(sat, getSatBCDESolar(sat)),
#'                            rada = getSatRADA(sat, getSatBCDESolar(sat)),
#'                            szen = getSatSZEN(sat, getSatBCDESolar(sat)),
#'                            esun = getSatESUN(sat, getSatBCDESolar(sat)),
#'                            model = "DOS2")
#' 
#' sensor_rad <- convSC2Rad(x = getSatDataLayer(sat, bcde), 
#'                          mult = getSatRADM(sat, bcde), 
#'                          add = getSatRADA(sat, bcde), getSatSZEN(sat, bcde))
#' 
#' ref_atmos <- calcAtmosCorr(x = sensor_rad,
#'                            path_rad = path_rad[names(path_rad) == bcde],
#'                            esun = getSatESUN(sat, bcde),
#'                            szen = getSatSZEN(sat, bcde), 
#'                            model = "DOS2")
#'                            
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return Satellite object with added atmospheric corrected layers
#' 
#' @rdname calcAtmosCorr
#'
setMethod("calcAtmosCorr", 
          signature(x = "Satellite"), 
          function(x, model = c("DOS2", "DOS4"), esun_method = "RadRef"){

            # If not supplied, 'model' defaults to DOS2
            model <- model[1]
            
            # Get solar bands with calibration information equals scaled counts
            sc_bands <- getSatBCDESolarCalib(x, calib = "SC")
            
            # Path radiance computation if necessary (TOA irradiance will be
            # computed, too, if necessary)
            if(any(is.na(getSatPRAD(x, sc_bands)))){
              x <- calcPathRadDOS(x, model = model, esun_method = esun_method)
            }
            
            # Radiance conversion if necessary
            if(any(is.na(getSatBCDESolarCalib(x, calib = "RAD")))){
              x <- convSC2Rad(x, szen_correction = "TRUE")
            }
            
            # Compute atmospheric correction (reflectance)
            rad_bands <- getSatBCDESolarCalib(x, calib = "RAD")
            for(bcde_rad in rad_bands){
              ref <-calcAtmosCorr(getSatDataLayer(x, bcde_rad),
                                  path_rad = getSatPRAD(x, bcde_rad),
                                  esun = getSatESUN(x, bcde_rad),
                                  szen = getSatSZEN(x, bcde_rad), 
                                  model = "DOS2")
              layer_bcde <- paste0(substr(bcde_rad, 1, nchar(bcde_rad) - 4),
                                   "_REF_AtmosCorr")
              meta_param <- data.frame(getSatSensorInfo(x),
                                       getSatBandInfo(x, bcde_rad, 
                                                      return_calib = FALSE),
                                       CALIB = "REF_AtmosCorr")
              info <- sys.calls()[[1]]
              info <- paste0("Add layer from ", info[1], "(", 
                             toString(info[2:length(info)]), ")")
              x <- addSatDataLayer(x, bcde = layer_bcde, data = ref,
                                   meta_param = meta_param,
                                   info = info, in_bcde = bcde_rad)
            }
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @return raster::RasterStack object with atmospheric corrected layers
#' 
#' @rdname calcAtmosCorr
#'
setMethod("calcAtmosCorr", 
          signature(x = "RasterStack"), 
          function(x, path_rad, esun, szen, model = c("DOS2", "DOS4")){
            # If not supplied, 'model' defaults to DOS2
            model <- model[1]
            
            for(l in seq(nlayers(x))){
              x[[l]] <- calcAtmosCorr(x, path_rad, esun, szen, model = "DOS2")
            }
            return(x)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @return raster::RasterLayer object with atmospheric corrected layer
#' 
#' @rdname calcAtmosCorr
#'
setMethod("calcAtmosCorr", 
          signature(x = "RasterLayer"), 
          function(x, path_rad, esun, szen, model = c("DOS2", "DOS4")){
            # If not supplied, 'model' defaults to DOS2
            model <- model[1]
            
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
            
            ref <- pi * (x - path_rad) / 
              ( tv * (esun * cos_szen * tz + edown) )
            
            return(ref)
          })          