#' Simulated Digital PCR data
#' 
#' Simulated data from array-based digital PCR experiment (see \code{\link{sim_adpcr}}).
#' 
#' @name six_panels
#' @docType data
#' @format An object of class \code{\linkS4class{adpcr}} containing six runs from three 
#' experiments (two runs per each experiment).
#' @examples 
#' #code below was used to create six_panels data set
#' \dontrun{
#' set.seed(1944)
#' adpcr1 <- sim_adpcr(m = 10, n = 765, times = 10000, pos_sums = FALSE, n_panels = 2)
#' adpcr2 <- sim_adpcr(m = 40, n = 765, times = 10000, pos_sums = FALSE, n_panels = 2)
#' adpcr2 <- rename_dpcr(adpcr2, exper = "Experiment2")
#' adpcr3 <- sim_adpcr(m = 100, n = 765, times = 10000, pos_sums = FALSE, n_panels = 2)
#' adpcr3 <- rename_dpcr(adpcr3, exper = "Experiment3")
#' bind_dpcr(adpcr1, adpcr2, adpcr3)
#' }
NULL