# Created by Matthew A. Birk
# Dependencies: None
# Converts txt file output from PreSens fiber optic O2 transmitters into usable dataframe
# Last updated: Mar 2016

#' Import data from PreSens O2 transmitter
#'
#' Imports the standard txt file output from most PreSens fiber optic O2 transmitters and converts the data into a data frame.
#'
#' The following PreSens fiber optic O2 transmitters are supported:
#' \describe{
#' \item{Fibox 3}{}
#' \item{Fibox 3 trace}{}
#' \item{Fibox 3 LCD trace}{}
#' \item{Microx TX3}{}
#' \item{Microx TX3 trace}{}
#' \item{OXY-4 mini}{}
#' \item{OXY-4 micro}{}
#' \item{OXY-4 trace}{}
#' \item{OXY-10 mini}{}
#' \item{OXY-10 micro}{}
#' \item{OXY-10 trace}{}
#' }
#' It is very important to note that the PreSens fiber optics O2 transmitters that are supported with this function DO NOT account for salinity (i.e. they assume salinity = 0 ppt). If the water sample measured was not fresh water, the oxygen concentrations (e.g. mg per liter or umol per liter) are incorrect in the PreSens txt file. This function corrects these O2 concentrations based on the salinity value defined by the \code{salinity} argument. Absolute partial pressures (i.e. hPa and torr) will also be slightly different due to the slight influence of salinity on water's vapor pressure. This difference is typically ~0.05\% of the recorded value.
#'
#' @param file a character string. The filepath for the file to be read.
#' @param o2_unit a character string. The unit of O2 measurement to be output in the data.frame. Options are:\describe{
#' \item{percent_a.s. (percent air saturation)}{}
#' \item{percent_o2}{}
#' \item{hPa}{}
#' \item{kPa}{}
#' \item{torr}{}
#' \item{mmHg}{}
#' \item{inHg}{}
#' \item{mg_per_l}{}
#' \item{umol_per_l}{}
#' \item{ml_per_l}{}
#' }
#' @param date a character string. The date format to be passed to \code{\link{strptime}}.
#' @param salinity salinity of water sample (psu). Default is 35 psu.
#'
#' @return A data frame with seven columns is returned.
#' \describe{
#' \item{TIME}{Date and time, POSIXct format.}
#' \item{DURATION}{Duration of measurement trial (minutes).}
#' \item{oxygen}{Oxygen measurement in desired unit. Column name changes based on \code{o2_unit} argument.}
#' \item{PHASE}{Phase recorded. Phase is inversely related to O2.}
#' \item{AMPLITUDE}{Amplitude recorded. Amplitude is an indicator of the quality of the signal. A low amplitude warning is produced by the transmitter below 2500.}
#' \item{TEMPERATURE}{Temperature recorded or defined at beginning of measurement trial.}
#' \item{ERROR_CODE}{Error code from transmitter. See PreSens user manual for translation of error code.}
#' }
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#'
#' @note Conversions are estimates based on the \code{\link[marelac]{marelac}} package and therefore differ slightly from the conversions provided by PreSens.
#' 
#' @examples
#' \dontrun{
#' file <- system.file('extdata', 'all_o2_units.txt', package = 'presens')
#' import_o2(file, o2_unit = 'umol_per_l', salinity = 25)
#' }
#'
#' @encoding UTF-8
#' @export
#' @import stats

import_o2 = function(file, o2_unit = 'percent_a.s.', date = '%d/%m/%y', salinity = 35)
{
  raw = readLines(file)
  raw = gsub(pattern = '\xb0|\xa9', replacement = ' ', raw) # replace non ASCII characters
  raw = raw[sapply(raw, nchar) > 0] # remove blank rows
  f = suppressWarnings(raw[grep('logtime|Logtime', raw):length(raw)]) # start dataframe from one row below the word 'logtime'
  f = gsub(pattern = ' ', replacement = '', f)
  f = strsplit(f, split = ';')
  f = as.data.frame(matrix(unlist(f), ncol = 8, byrow = TRUE), stringsAsFactors = FALSE)
  o2_col = ifelse(o2_unit %in% names(o2_unit_conv()), toupper(o2_unit), stop('the o2_unit argument is not an acceptable unit', call. = F))
  air_pres = suppressWarnings(as.numeric(stats::na.omit(as.numeric(unlist(strsplit(raw[grep('Pressure', raw)], ' |;')))))/1000)
  colnames(f) = c('DATE', 'TIME', 'DURATION', o2_col, 'PHASE', 'AMPLITUDE', 'TEMPERATURE', 'ERROR_CODE')
  unit_id_index = grep('oxygen', f[, o2_col])
  o2_string_options = list(
    percent_a.s. = 'oxygen/%airsatur.',
    percent_o2 = 'oxygen/%O2',
    hPa = 'oxygen/hPa(mbar)',
    torr = 'oxygen/Torr',
    mg_per_l = 'oxygen/mg/L(ppm)',
    umol_per_l = 'oxygen/umol/L'
  )
  f_split = split(f, findInterval(1:nrow(f), unit_id_index))
  f_split = lapply(f_split, function(i){
  	inter1 = names(which(i[1, o2_col] == o2_string_options))
  	i = i[-1, ]
  	i[, o2_col] = as.numeric(i[, o2_col])
  	i[, o2_col] = i[, o2_col] / o2_unit_conv(to = inter1, salinity = 0, temp = as.numeric(i$TEMPERATURE), air_pres = air_pres)[[inter1]] * o2_unit_conv(to = o2_unit, salinity = salinity, temp = as.numeric(i$TEMPERATURE), air_pres = air_pres)[[o2_unit]]
  	i
  })
	f = do.call('rbind', f_split)
  f[, 'TIME'] = as.data.frame(strptime(paste(f$DATE, f$TIME), paste(date, '%T')))
  f$DATE = NULL
  f$DURATION = as.numeric(f$DURATION)
  if(any(is.na(f$TIME))) stop(paste('The time record does not match', date, 'on at least some of the lines between', range(which(is.na(f$TIME)))[1], 'and', range(which(is.na(f$TIME)))[2]), call. = F)
  f$PHASE = as.numeric(f$PHASE)
  f$AMPLITUDE = as.numeric(f$AMPLITUDE)
  f$TEMPERATURE = as.numeric(f$TEMPERATURE)
  f$ERROR_CODE = as.factor(f$ERROR_CODE)
  if(length(unit_id_index) > 1) message('Multiple units of O2 were recorded. All values were converted to ', o2_unit, call. = F)
  if(any(f$ERROR_CODE != 'E0')) warning('Errors recorded during trial! Check ERROR_CODE column', call. = F)
  row.names(f) = NULL
  return(f)
}