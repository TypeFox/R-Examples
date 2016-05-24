#' @name watbal.simobsdata
#' @title Soil water content measurements and associated simulations with WaterBalance model
#' @description
#' Data of soil water content from Luc Champolivier (CETIOM), in En Crambade (31, France), on canola without irrigation in 2008.
#' sonde Diviner 2000 (from Sentek Pty Ltd)
#' Simulation are from watbal.model, with an initial water content estimated from measurement with Diviner 2000.
#' @docType data
#' @usage watbal.simobsdata
#' @format a \code{RangedData} instance, 1 row per day :
#' Weather : day / RAIN / ETr /
#' simulation :  WAT / WATp / ARID
#' observation :t1_WATp_0_40cm / t2_WATp_0_40cm / t3_WATp_0_40cm / WATp_SF.mean / WATp_SF.var
#' @source CETIOM, Luc Champolivier (2008).
NULL
