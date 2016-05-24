#' Get calibration information from Landsat 8 standard level 1B/T filename
#'
#' @description
#' The function scans a Lansat metadata file for various calibration 
#' and orbit coefficients as well as some sensor specific data.
#'
#' @param files Path and filename of the Landsat metadata file. 
#'
#' @return \code{data.frame} containing the following information for each 
#' band/layer: 
#' \itemize{
#'   \item DATE date (e.g. 2013-07-07)
#'   \item SID sensor id (e.g. LC8)
#'   \item SENSOR sensor name (e.g. Landsat 8)
#'   \item SGRP sensor group (e.g. Landast)
#'   \item BID band id (e.g. 7)
#'   \item BCDE band code (5 digit standard name, e.g B001n)
#'   \item SRES spatial resolution of the sensor band (e.g. 30 for 30 m x 30m)
#'   \item TYPE type of the sensor band regarding wavelength (e.g. VIS)
#'   \item SPECTRUM spectral range regarding radiation source (e.g. solar)
#'   \item CALIB type of applied calibration (e.g. SC for scaled counts)
#'   \item RADA addtition coefficient for radiance conversion
#'   \item RADM multiplication coefficient for radiance conversion
#'   \item REFA addtition coefficient for reflectance conversion
#'   \item REFM multiplication coefficient for reflectance conversion
#'   \item BTK1 brightness temperature correction parameter
#'   \item BTK2 brightness temperature correction parameter
#'   \item SZEN sun zenith angle
#'   \item SAZM sun azimuth angle
#'   \item SELV sun elevation angle
#'   \item ESD earth-sun distance (AU)
#'   \item LMIN Minimum wavelength of the band (micrometer)
#'   \item LMAX Maximum wavelength of the band (micrometer)
#'   \item RADMIN Minimum radiance recorded by the band
#'   \item RADMAX Maximum radiance recorded by the band
#'   \item REFMIN Minimum reflectance recorded by the band
#'   \item REFMAX Maximum reflectance recorded by the band
#'   \item LNBR Layer number from 1 to n layers
#'   \item LAYER Layer name
#'   \item FILE Filepath of the data file
#'   \item METAFILE Filepath of the metadata file
#' }
#' 
#' @export compMetaLandsat
#'
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' compMetaLandsat(files)
#' 
compMetaLandsat <- function(files){
  
  datafiles <- compFilePathLandsat(files)
  
  bandinfo <- lutInfoBandsFromSID(datafiles$SID[1])
  bandinfo <- merge(bandinfo, datafiles, by = "BCDE")
  
  metadata <- utils::read.table(as.character(bandinfo$METAFILE[1]), header = FALSE, 
                         sep = "=", fill = TRUE)
  
  search_term_date <- "DATE_ACQUIRED"
  search_term_esd <- "EARTH_SUN_DISTANCE"
  
  metainformation <- lapply(seq(nrow(bandinfo)), function(x){
    search_term_rad_max <- paste0("RADIANCE_MAXIMUM_BAND_", x)
    search_term_rad_min <- paste0("RADIANCE_MINIMUM_BAND_", x)
    search_term_ref_max <- paste0("REFLECTANCE_MAXIMUM_BAND_", x)
    search_term_ref_min <- paste0("REFLECTANCE_MINIMUM_BAND_", x)
    
    search_term_add_rad <- paste0("RADIANCE_ADD_BAND_", x)
    search_term_mult_rad <- paste0("RADIANCE_MULT_BAND_", x)
    search_term_add_ref <- paste0("REFLECTANCE_ADD_BAND_", x)
    search_term_mult_ref <- paste0("REFLECTANCE_MULT_BAND_", x)
    search_term_BTK1 <- paste0("K1_CONSTANT_BAND_", x)
    search_term_BTK2 <- paste0("K2_CONSTANT_BAND_", x)
    
    cal_rad_max <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_rad_max))))
    cal_rad_min <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_rad_min))))
    cal_ref_max <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_ref_max))))
    cal_ref_min <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_ref_min))))
    
    cal_add_rad <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_add_rad))))
    cal_mult_rad <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_mult_rad))))
    cal_add_ref <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_add_ref))))
    cal_mult_ref <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_mult_ref))))
    cal_BTK1 <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_BTK1))))
    cal_BTK2 <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_BTK2))))
    date <- as.POSIXlt(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_date))), 
      tz = "UTC")
    esd <- as.numeric(as.character(
      (subset(metadata$V2, gsub("\\s","", metadata$V1) == search_term_esd))))
    selv <- as.numeric(as.character(
      subset(metadata$V2, gsub("\\s","", metadata$V1) == "SUN_ELEVATION")))
    sazm <- as.numeric(as.character(
      subset(metadata$V2, gsub("\\s","", metadata$V1) == "SUN_AZIMUTH")))
    szen <- 90.0 - selv

    if(length(cal_rad_max) == 0){cal_rad_max = NA}
    if(length(cal_rad_min) == 0){cal_rad_min = NA}
    if(length(cal_ref_max) == 0){cal_ref_max = NA}
    if(length(cal_ref_min) == 0){cal_ref_min = NA}
    if(length(cal_add_rad) == 0){cal_add_rad = NA}
    if(length(cal_mult_rad) == 0){cal_mult_rad = NA}
    if(length(cal_add_ref) == 0){cal_add_ref = NA}
    if(length(cal_mult_ref) == 0){cal_mult_ref = NA}
    if(length(cal_BTK1) == 0){cal_BTK1 = NA}
    if(length(cal_BTK2) == 0){cal_BTK2 = NA}
    if(length(esd) == 0){esd = NA}
    
    result <- data.frame(DATE = date, 
                         SID = bandinfo$SID.x[x],
                         SENSOR = bandinfo$SENSOR[x],
                         SGRP = bandinfo$SGRP[x],
                         BID = bandinfo$BID.x[x],
                         BCDE = bandinfo$BCDE[x],
                         #SRES = bandinfo$SRES[x],
                         TYPE = bandinfo$TYPE[x],
                         SPECTRUM = bandinfo$SPECTRUM[x],
                         CALIB = bandinfo$CALIB[x],
                         RADA = cal_add_rad,
                         RADM = cal_mult_rad,
                         REFA = cal_add_ref,
                         REFM = cal_mult_ref,
                         BTK1 = cal_BTK1,
                         BTK2 = cal_BTK2,
                         SZEN = szen,
                         SAZM = sazm,
                         SELV = selv,
                         ESD = esd,
                         LMIN = bandinfo$LMIN[x],
                         LMAX = bandinfo$LMAX[x],
                         RADMAX = cal_rad_max,
                         RADMIN = cal_rad_min,
                         REFMAX = cal_ref_max,
                         REFMIN = cal_ref_min,
                         LNBR = x,
                         LAYER = bandinfo$LAYER[x],
                         FILE = bandinfo$FILE[x],
                         METAFILE = bandinfo$METAFILE[x],
                         stringsAsFactors = FALSE)
  })
  metainformation <- do.call("rbind", metainformation)
  return(metainformation)
}