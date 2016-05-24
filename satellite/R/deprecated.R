#' Deprecated functions
#' 
#' @name deprecated
#' 
#' @param x Satellite object
#' @param convert Convert dataset to radiance/reflectance/temperature
#' @param szen_correction Compute sun zenith correction
#' @param band Band of the sensor
#' @param mult Multiplicative calibration coefficient
#' @param add Additive multiplication coefficient
#' @param szen Sun zenith angle
#' @param k1 BTT calibration coefficient
#' @param k2 BTT calibration coefficient
#' @param method Method used for computation
#' @param model Model used for computation 
#' @param normalize Normalize earth-sun distance
#' @param esd Earth sun distance
#' @param atmos_model Atmospheric model to be used
#' @param esun_mode Earth sun distance computation method to be used
#' @param ref_mult Multiplicative reflection calibration coefficient
#' @param ref_add Additive reflection calibration coefficient
#' @param rad_mult Multiplicative radiance calibration coefficient
#' @param rad_add Additive radiance calibration coefficient
#'
#' @description
#' The functions have been implemented in the very beginning of the package
#' development, mainly to be used within the scope of a remote sensing course at 
#' Marburg University. To ensure that the scripts developed within this course
#' will still work after the next major revision, they are still part of this
#' package, but they will mainly just foreward the respective call to the
#' more up-to-date function.
#' 
NULL

# Deprecated satCalib ----------------------------------------------------------
#' @export satCalib
#'
#' @rdname deprecated
#'
satCalib <- function(x, convert = "all", szen_correction = "TRUE"){
  .Deprecated("convDN2RU")
  convDN2RU(x, convert, szen_correction)
}


# Deprecated calibLinear -------------------------------------------------------
#' @export calibLinear
#'
#' @rdname deprecated
#'
calibLinear <- function(band, mult, add, szen, k1, k2){
  .Deprecated("convDN2RU")
  if(missing(szen)){
    if(missing(k1)){
      convDN2RU(x = band, mult = mult, add = add)
    } else {
      convDN2RU(x = band, mult = mult, add = add, k1 = k1, k2 =k2)  
    }
  } else {
    if(missing(k1)){
      convDN2RU(x = band, mult = mult, add = add, szen = szen)
    } else {
      convDN2RU(x = band, mult = mult, add = add, szen = szen,
                      k1 = k1, k2 =k2)
    }
  }
}


# Deprecated satTOAIrrad -------------------------------------------------------
#' @export satTOAIrrad
#'
#' @rdname deprecated
#'
satTOAIrrad <- function(x, method = "Table", model = "MNewKur", 
                        normalize = TRUE, esd){
  if((method != "RadRef" & normalize == FALSE & missing(esd)) |
       (method == "RadRef" & normalize == TRUE & missing(esd))){
    esd = getSatESD(x)
    if(is.na(esd)){
      esd = calcEarthSunDist(date)
    } 
  }
  if(method == "Table"){
    if(normalize == TRUE){
      esun <- calcTOAIrradTable(x = getSatSID(x), 
                                normalize  = normalize)
    } else {
      esun <- calcTOAIrradTable(x = getSatSID(x), 
                                normalize  = normalize, 
                                esd = esd)
    }
    bcde = names(esun)
  } else if(method == "Model"){
    rsr <- lutInfoRSRromSID(sid = getSatSID(x))
    if(normalize == TRUE){
      esun <- calcTOAIrradModel(x = rsr, model = model, 
                                normalize = normalize)
    } else {
      esun <- calcTOAIrradModel(x = rsr, model = model, 
                                normalize = normalize, esd = esd)
    }
    bcde = names(esun)
  } else if(method == "RadRef"){
    if(normalize == TRUE){
      esun <- 
        calcTOAIrradRadRef(
          x = getSatRadMax(x, getSatBCDESolar(x)), 
          ref_max = getSatRefMax(x, getSatBCDESolar(x)),
          esd = esd, normalize = normalize)
    } else {
      esun <- 
        calcTOAIrradRadRef(
          x = getSatRadMax(x, getSatBCDESolar(x)), 
          ref_max = getSatRefMax(x, getSatBCDESolar(x)), 
          normalize = normalize)
    }
    bcde = getSatBCDESolar(x)
  }
  x <- addSatMetaParam(x, meta_param = data.frame(BCDE = bcde,
                                                  ESUN = as.numeric(esun)))
  return(x)
}


# Deprecated satTOAIrrad -------------------------------------------------------
#' @export satPathRadDOS
#'
#' @rdname deprecated
#'
satPathRadDOS <- function(x, atmos_model = "DOS2", esun_mode = "RadRef"){
  calcPathRadDOS(x, model = atmos_model, esun_method = esun_mode)
}


# Deprecated satAtmosCorr ------------------------------------------------------
#' @export satAtmosCorr
#'
#' @rdname deprecated
#'
satAtmosCorr <- function(x, atmos_model = "DOS2", esun_mode = "RadRef"){
  calcAtmosCorr(x, model = atmos_model, esun_method = esun_mode)
}


# Deprecated calibLinearInverse ------------------------------------------------
#' @export calibLinearInverse
#'
#' @rdname deprecated
#'
calibLinearInverse <- function(band, ref_mult, ref_add, rad_mult, rad_add, 
                               szen){
  if(missing(szen)){
    convRef2RadLinear(band, ref_mult, ref_add, rad_mult, rad_add)
  } else {
    convRef2RadLinear(band, ref_mult, ref_add, rad_mult, rad_add, szen)  
  }
}


# Deprecated satInvarFeatures --------------------------------------------------
#' @export satInvarFeatures
#'
#' @rdname deprecated
#'
satInvarFeatures <- function(x){
  maskInvarFeatures(x)
}