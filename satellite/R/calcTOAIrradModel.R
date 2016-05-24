if ( !isGeneric("calcTOAIrradModel") ) {
  setGeneric("calcTOAIrradModel", function(x, ...)
    standardGeneric("calcTOAIrradModel"))
}
#' Compute top of atmosphere solar irradiance for sensor bands using LUTs
#'
#' @description
#' Compute mean extraterrestrial solar irradiance (ESun) using tabulated mean
#' solar spectral data and the band specific relative spectral response (rsr) 
#' functions.
#' 
#' @param x A Satellite object or the relative spectral response function for 
#' the respective band as \code{data.frame} (see details for structure).
#' @param model Tabulated solar radiation model to be used (one of 
#' MCebKur_MChKur, MNewKur,  MthKur, MoldKur, MODWherli_WMO, NN, see reference
#' on tabulated solar irradiance below).
#' @param normalize Logical; if \code{TRUE}, ESun is normalized to mean 
#' earth-sun distance. 
#' @param esd Earth-sun distance (AU, can be estimated using 
#' \code{\link{calcEarthSunDist}}). If x is a Satellite object and esd is not 
#' supplied and necessary for normalization, it is tried to take it from the 
#' metadata, otherwise it is estimated by the day of the year using 
#' \code{\link{calcEarthSunDist}}.
#'
#' @export calcTOAIrradModel
#' 
#' @name calcTOAIrradModel
#' 
#' @details Computation of ESun is taken from Updike and Comp (2011).
#' 
#' Tabulated values for mean earth-sun distance are taken from the 
#' data sources mentioned in the references. 
#' 
#' If results should not be normalized to a mean earth-sun distance, the 
#' actual earth-sun distance is approximated by the day of the year using
#' \code{\link{calcEarthSunDist}}.
#' 
#' Relative spectral response values have to be supplied as a \code{data.frame}
#' which has at least the following three columns: (i) a column "Band" for the 
#' sensor band number (i.e. 1, 2, etc.), (ii) a column "WAVELENGTH" for the 
#' WAVELENGTH data in full nm steps, and (iii) a column "RSR" for the response 
#' information [0...1].
#' 
#' @references 
#' Updike T, Comp C (2011) Radiometric use of WorldView-2 imagery. 
#' Technical Note, available online at \url{https://goo.gl/np4Bwm}.
#' 
#' Tabulated relative spectral response functions (nm-1) are taken from taken from the
#' \href{http://landsat.usgs.gov/instructions.php}{spectral viewer}
#' of the USGS Landsat FAQ.
#' 
#' Tabulated solar irradiance (W m-2 nm-1) is taken from the 
#' \href{http://rredc.nrel.gov/solar/spectra/am0/modtran.html}{National Renewable 
#' Energy Laboratory}.
#' 
#' @seealso \code{\link{calcTOAIrradTable}} for tabulated solar irradiance
#' values from the literature or \code{\link{calcTOAIrradRadRef}} for the 
#' computation of the solar irradiance based on maximum radiation and reflection
#' values of the dataset.
#' 
#' See \code{\link{calcEarthSunDist}} for calculating the earth-sun
#' distance based on the day of the year which is called by this function if
#' ESun should be corrected for actual earth-sun distance.
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' sat <- calcTOAIrradModel(sat)
#' getSatESUN(sat)
#' 
#' lut <- lutInfo()
#' calcTOAIrradModel(lut$L8_RSR, model = "MNewKur", normalize = FALSE, 
#'   esd = calcEarthSunDist("2015-01-01"))
#' 
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return If x is a satellite object, a Satellite object with ESun information 
#' added to the metadata; if x is a \code{data.frame}, a vector containing ESun
#' for the respective band(s).
#' 
#' @rdname calcTOAIrradModel
#'
setMethod("calcTOAIrradModel", 
          signature(x = "Satellite"), 
          function(x, model = "MNewKur", normalize = TRUE, esd){
            if(normalize == FALSE & missing(esd)){
              esd <- getSatESD(x)
              if(is.na(esd)){
                doy <- getSatDATE(x)
                esd <- calcEarthSunDist(doy)
              } 
            }
            rsr <- lutInfoRSRromSID(sid = getSatSID(x))
            if(normalize == TRUE){
              esun <- calcTOAIrradModel(x = rsr, model = model, 
                                        normalize = normalize)
            } else {
              esun <- calcTOAIrradModel(x = rsr, model = model, 
                                        normalize = normalize, esd = esd)
            }
            bcde <- names(esun)
            
            x <- addSatMetaParam(x, 
                                 meta_param = data.frame(
                                   BCDE = bcde,
                                   ESUN = as.numeric(esun)))
            return(x)
          })


# Function using data frame ----------------------------------------------------
#' 
#' @rdname calcTOAIrradModel
#'
setMethod("calcTOAIrradModel", 
          signature(x = "data.frame"), 
          function(x, model = "MNewKur", normalize = TRUE, esd){
            toa <- lut$SOLAR
            toa$WAVELENGTH <- round(toa$WAVELENGTH, 0)
            toa_aggregated <- aggregate(toa, by = list(toa$WAVELENGTH), 
                                        FUN = "mean")
            
            rsr <- x[, c(grep("BCDE", colnames(x)), 
                           grep("WAVELENGTH", colnames(x)), 
                           grep("RSR", colnames(x)))]
            
            eSun <- lapply(unique(rsr$BCDE), function(y){
              rsr_solar <- merge(rsr[rsr$BCDE == y,], toa_aggregated, 
                                 by = "WAVELENGTH")
              if(nrow(rsr_solar) > 0){
                act_eSun <- aggregate(
                  rsr_solar$RSR * rsr_solar[,grep(model, names(rsr_solar))], 
                  by = list(rsr_solar$BCDE), 
                  FUN = "sum")[2] / aggregate(rsr_solar$RSR, 
                                              by = list(rsr_solar$BCDE), 
                                              FUN = "sum")[2] * 1000
                act_eSun <- unlist(act_eSun)
              } else {
                act_eSun <- c(NA)
              }
              names(act_eSun) <- paste0(y)
              return(act_eSun)
            })
            eSun <- unlist(eSun)
            
            if(normalize == FALSE){
              if(missing(esd)){
                stop("Variable esd is missing.")
              }
              eSun <- eSun / esd**2
            }
            return(eSun)
          })