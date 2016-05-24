#' Find an optimal p-value for SONE
#'
#' @description a wrapper that runs the maximum and minimum scenarios using \code{\link{scenario_sim}} and provides the optimal p -value
#' @inheritParams scenario_sim
#' @param n_indicators How many many indicators are there in a scale. The package is tested with 8 indicators (default), but should work with other number.
#' @param plotting Plots the result with \code{\link{optimal_p_out}}. Defaults to ''. Possible options: '' - no plot; 'yes' - a regular plot;  'file' -- writes the plot to a tiff file in working directory. If sizes is a single value, plotting is disabled.
#' @param to_min How many indicators relate to the outcome in the lack of ION condition. In \code{\link{optimal_p}} defaults to (round((n_indicators/2),0)) - 1), i.e close to half the number of indicators.

#' @return Returns the P criterion, as well as the p values for max and min scenario for each sample size. If min pvalue > max pvalue, then p criterion is NA.
#' @encoding utf-8

#' @examples
#' set.seed(466)
#' n_sim=100
#' ptm <- proc.time()
#' a=optimal_p(sizes=750, n_sim=n_sim, n_indicators=8, cor_to_outcome=0.25)
#' stp=proc.time() - ptm
#' print(paste("Currently elapsed:",round(stp[3],1)))
#' print(paste("Time estimate for n_sim=5000:",round(stp[3]*5000/n_sim,1)))
#'
#' @export



optimal_p <- function(sizes, n_sim = 100, plotting = "", n_indicators = 8, to_min = (round((n_indicators/2), 0)) -
    1, ...) {

    scenario_max <- scenario_sim(sizes = sizes, n_sim = n_sim, to_n = n_indicators, tn_n = 0, ...)
    scenario_min <- scenario_sim(sizes = sizes, n_sim = n_sim, to_n = to_min, tn_n = n_indicators - to_min, ...)
    tab = optimal_p_out(scenario_max = scenario_max[[1]], scenario_min = scenario_min[[1]], to_min = to_min, sizes = sizes,
        n_sim = n_sim, plotting = plotting)
    out = list(tab, scenario_max, scenario_min)
    return(out)
}

