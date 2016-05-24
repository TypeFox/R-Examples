# Created by Matthew A. Birk
# Dependencies: birk, marelac
# Converts % air saturation to other O2 units
# Last updated: Mar 2016

#' Convert units of dissolved oxygen
#'
#' Given a measurement of dissolved O2, a list of commonly used units of oxygen partial pressures and concentrations are returned.
#'
#' Conversions are based on relationships and values from the package \code{\link[marelac]{marelac}}.
#'
#' @param o2 a numeric vector of the O2 value(s). Default is 100.
#' @param from a string describing the unit used to measure \code{o2}. Default is "percent_a.s." Options are:\describe{
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
#' @param to a single string either describing the unit to which the conversion should be conducted (options are the same as in \code{from}), or the string "all" to return all units.
#' @param salinity salinity of water sample (psu). Default is 35 psu.
#' @param temp temperature of water sample (°C). Default is 25 °C.
#' @param air_pres pressure of air overlying water sample (bar). Default is 1.013253 bar.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#'
#' @examples
#' o2_unit_conv(o2 = 50)
#' o2_unit_conv(o2 = 1:50, from = "umol_per_l", to = "ml_per_l", salinity = 0, temp = 10,
#' 	air_pres = 1.2)
#' o2_unit_conv()[c('mmHg','kPa')]
#'
#' @encoding UTF-8
#' @export
#' @import marelac
#' @import birk

o2_unit_conv = function(o2 = 100, from = 'percent_a.s.', to = 'all', salinity = 35, temp = 25, air_pres = 1.013253){
	if(from == 'percent_a.s.') perc_a.s. = o2
	if(from == 'percent_o2') perc_a.s. = o2 / marelac::atmComp('O2')
	if(from == 'hPa') perc_a.s. = birk::conv_unit(o2, 'hPa', 'atm') * 100 / (air_pres - marelac::vapor(S = salinity, t = temp)) / marelac::atmComp('O2')
	if(from == 'kPa') perc_a.s. = birk::conv_unit(o2, 'kPa', 'atm') * 100 / (air_pres - marelac::vapor(S = salinity, t = temp)) / marelac::atmComp('O2')
	if(from == 'torr') perc_a.s. = birk::conv_unit(o2, 'torr', 'atm') * 100 / (air_pres - marelac::vapor(S = salinity, t = temp)) / marelac::atmComp('O2')
	if(from == 'mmHg') perc_a.s. = birk::conv_unit(o2, 'mmHg', 'atm') * 100 / (air_pres - marelac::vapor(S = salinity, t = temp)) / marelac::atmComp('O2')
	if(from == 'inHg') perc_a.s. = birk::conv_unit(o2, 'inHg', 'atm') * 100 / (air_pres - marelac::vapor(S = salinity, t = temp)) / marelac::atmComp('O2')
	if(from == 'mg_per_l') perc_a.s. = o2 * 100 / marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2') / 1e-6 / marelac::molweight('O2') / 1e3
	if(from == 'umol_per_l') perc_a.s. = o2 * 100 / marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2')
	if(from == 'ml_per_l') perc_a.s. = birk::conv_unit(o2, 'ml', 'l') / marelac::molvol(t = temp, P = air_pres, species = 'O2', quantity = 1 / birk::conv_unit(100 / marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2'), 'mol', 'umol'))
	
	x=list()
	if(to == 'percent_a.s.' | to == 'all') x$percent_a.s. = perc_a.s.
	if(to == 'percent_o2' | to == 'all') x$percent_o2 = marelac::atmComp('O2') * perc_a.s.
	if(to == 'hPa' | to == 'all') x$hPa = birk::conv_unit((air_pres - marelac::vapor(S = salinity, t = temp)) * marelac::atmComp('O2') * perc_a.s. / 100, 'atm', 'hPa')
	if(to == 'kPa' | to == 'all') x$kPa = birk::conv_unit((air_pres - marelac::vapor(S = salinity, t = temp)) * marelac::atmComp('O2') * perc_a.s. / 100, 'atm', 'kPa')
	if(to == 'torr' | to == 'all') x$torr = birk::conv_unit((air_pres - marelac::vapor(S = salinity, t = temp)) * marelac::atmComp('O2') * perc_a.s. / 100, 'atm', 'torr')
	if(to == 'mmHg' | to == 'all') x$mmHg = birk::conv_unit((air_pres - marelac::vapor(S = salinity, t = temp)) * marelac::atmComp('O2') * perc_a.s. / 100, 'atm', 'mmHg')
	if(to == 'inHg' | to == 'all') x$inHg = birk::conv_unit((air_pres - marelac::vapor(S = salinity, t = temp)) * marelac::atmComp('O2') * perc_a.s. / 100, 'atm', 'inHg')
	if(to == 'mg_per_l' | to == 'all') x$mg_per_l = marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2') * 1e-6 * marelac::molweight('O2') * 1e3 * perc_a.s. / 100
	if(to == 'umol_per_l' | to == 'all') x$umol_per_l = marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2') * perc_a.s. / 100
	if(to == 'ml_per_l' | to == 'all') x$ml_per_l = birk::conv_unit(as.numeric(marelac::molvol(t = temp, P = air_pres, species = 'O2', quantity = birk::conv_unit(marelac::gas_satconc(S = salinity, t = temp, P = air_pres, species = 'O2') * perc_a.s. / 100, 'umol', 'mol'))), 'l', 'ml')
	attr(x$percent_o2, 'names') = NULL
	attr(x$hPa, 'names') = NULL
	attr(x$kPa, 'names') = NULL
	attr(x$torr, 'names') = NULL
	attr(x$mmHg, 'names') = NULL
	attr(x$inHg, 'names') = NULL
	attr(x$mg_per_l, 'names') = NULL
	attr(x$umol_per_l, 'names') = NULL
	attr(x$ml_per_l, 'names') = NULL
	return(x)
}