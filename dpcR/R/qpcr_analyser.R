#' qPCR Analyser
#' 
#' Calculate statistics based on fluorescence. The function can be used to
#' analyze amplification curve data from quantitative real-time PCR
#' experiments. The analysis includes the fitting of the amplification curve by
#' a non-linear function and the calculation of a quantification point (often
#' referred to as Cp (crossing-point), Cq or Ct) based on a user defined
#' method. The function can be used to analyze data from chamber based dPCR
#' machines.
#' 
#' The \code{qpcRanalyzer} is a functions to automatize the analysis of
#' amplification curves from conventional quantitative real-time PCR (qPCR)
#' experiments and is adapted for the needs in dPCR.  This function calls
#' instances of the \code{qpcR} package to calculate the quantification
#' points (cpD1, cpD2, Cy0 (default), TOP (optional)), the amplification
#' efficiency, fluorescence at the quantification point (Cq), the absolute
#' change of fluorescence and the take-off point (TOP). Most of the central
#' functionality of the \code{qpcR} package is accessible. The user can assign
#' concentrations to the samples. One column contains binary converted (pos (1)
#' and neg (0)) results for the amplification reaction based on a user defined
#' criteria (Cq-range, fluorescence cut-off, ...). \code{qpcr_analyser} tries
#' to detect cases where an amplification did not take place of was impossible
#' to analyze. By default \code{qpcr_analyser} analyses uses the Cy0 as
#' described in Guescini et al. (2008) for estimation of the quantification
#' point since method is considered to be better suited for many probe systems.
#' By default a 5-parameter model is used to fit the amplification curves. As
#' such \code{qpcr_analyser} is a function, which serves for preliminary data
#' inspection (see Example section) and as input for other R functions from the
#' \code{dpcR} package (e.g., \link{plot_panel}).
#' 
#' @name qpcr_analyser
#' @aliases qpcr_analyser qpcr_analyser-methods qpcr_analyser,adpcr-method
#' qpcr_analyser,data.frame-method qpcr_analyser,modlist-method qpcr_analyser
#' @docType methods
#' @param input a dataframe containing the qPCR data or a result of function
#' \code{\link[qpcR]{modlist}} or an object of the class
#' \code{\linkS4class{adpcr}}.
#' @param cyc the column containing the cycle data. Defaults to first column.
#' @param fluo the column(s) (runs) to be analyzed. If NULL, all runs will be
#' considered.  Use fluo = 2 to chose the second column for example.
#' @param model is the model to be used for the analysis for all runs. Defaults
#' to 'l5' (see \code{\link[qpcR]{pcrfit}}).
#' @param norm logical. Indicates if the raw data should be normalized within
#' [0, 1] before model fitting.
#' @param iter_tr \code{iter_tr} number of iteration to fit the curve.
#' @param type is the method for the crossing point/threshold cycle estimation
#' and efficiency estimation (\link[qpcR]{efficiency}). Defaults to 'Cy0'
#' (\code{\link[qpcR]{Cy0}}).
#' @param takeoff logical; if \code{TRUE} calculates the first significant
#' cycle of the exponential region (takeoff point). See
#' \code{\link[qpcR]{takeoff}} for details.
#' @return A matrix where each column represents crossing point, efficiency,
#' the raw fluorescence value at the point defined by type and difference
#' between minimum and maximum of observed fluorescence. If takeoff parameter
#' is \code{TRUE}, additional two column represents start and the end of the
#' fluorescence growth.
#' @author Stefan Roediger, Andrej-Nikolai Spiess, Michal Burdukiewicz.
#' @seealso \link[qpcR]{modlist}.
#' @export
#' @references Ritz C, Spiess An-N, \emph{qpcR: an R package for sigmoidal
#' model selection in quantitative real-time polymerase chain reaction
#' analysis}.  Bioinformatics 24 (13), 2008.
#' 
#' Andrej-Nikolai Spiess (2013). qpcR: Modelling and analysis of real-time PCR
#' data.\cr \url{http://CRAN.R-project.org/package=qpcR}\cr
#' @keywords qPCR Cy0 real-time amplification quantification
#' @examples
#' 
#' # Take data of guescini1 data set from the qpcR R package.
#' library(qpcR)
#' # Use the first column containing the cycles and the second column for sample F1.1.
#' data(guescini1)
#' qpcr_analyser(guescini1, cyc = 1, fluo = 2)
#' 
#' # Use similar setting as before but set takeoff to true for an estimation of
#' # the first significant cycle of the exponential region.
#' qpcr_analyser(guescini1, cyc = 1, fluo = 2, takeoff = TRUE)
#' 
#' # Use similar setting as before but use qpcr_analyser in a loop to calculate the results for the
#' # first four columns containing the fluorescence in guescini1
#' print(qpcr_analyser(guescini1, cyc = 1, fluo = 2:5, takeoff = TRUE))
#' 
#' # Run qpcr_analyser on the list of models (finer control on fitting model process)
#' models <- modlist(guescini1)
#' qpcr_analyser(models)
#' 
NULL


qpcr_analyser <- function(input, cyc = 1, fluo = NULL, model = l5, norm = FALSE, iter_tr = 50, 
                           type = "Cy0", takeoff = FALSE) {
  stop("Wrong class of 'input'")
}

setGeneric("qpcr_analyser")

setMethod("qpcr_analyser", signature(input = "data.frame"), function(input, 
                                                                     cyc = 1, 
                                                                     fluo = NULL, 
                                                                     model = l5, 
                                                                     norm = FALSE, 
                                                                     iter_tr = 50, 
                                                                     type = "Cy0", 
                                                                     takeoff = FALSE) {
  all_fits <- fit_adpcr(input, cyc, fluo, model, norm, iter_tr)
  res <- analyze_qpcR(all_fits, type, takeoff)
  res <- cbind(res, deltaF = calc_deltaF(input, cyc, fluo))
  res
})


setMethod("qpcr_analyser", signature(input = "modlist"), function(input, type = "Cy0", takeoff = FALSE) {
  res <- analyze_qpcR(input, type, takeoff)
  res
})

analyze_qpcR <- function(fit_list, type = "Cy0",  takeoff = FALSE) {
  res <- t(vapply(fit_list, function(fit) 
    safe_efficiency(fit, type), c(0, 0, 0)))
  if (takeoff) {
    res <- cbind(res, t(vapply(fit_list, function(fit) 
      unlist(takeoff(fit)[c("top", "f.top")]), c(0, 0))))
  }
  
  res
}

calc_deltaF <- function(pcr_data, cyc, fluo) {
  if (is.null(fluo)) {
    all_fluos <- (1L:ncol(pcr_data))[-cyc]
  } else {
    all_fluos <- fluo
  }
  vapply(all_fluos, function(x)  
    quantile(tail(pcr_data[, x]), 0.85) - quantile(head(pcr_data[, x]), 0.25), 0)
}

#to do
#general function to test dpcr objects
#k - positive droplets/chambers
#total number of droplets/chambers
#different principle than dube/bhat method. Calculate confidence interval for k, not for m
# test_dpcr <- function(k, n) {
#   #theoretical values
#   theor <- round(fl(unlist(binom.confint(k, n, methods = "wilson", conf.level = 0.95)[, 4:6]))*n, 0)
#   sim_dpcr(theor[1], n, times = 1000, dube = TRUE, pos_sums = TRUE, n_panels = 1000)
# }


setMethod("qpcr_analyser", signature(input = "adpcr"), function(input, cyc = 1, fluo = NULL, 
                                                                model = l5, 
                                                                norm = FALSE, iter_tr = 50, 
                                                                type = "Cy0", takeoff = FALSE) {
  if (slot(input, "type") != "fluo")
    stop("'input' must contain fluorescence data.", call. = TRUE, domain = NA)
  input <- slot(input, ".Data")
  all_fits <- fit_adpcr(input, cyc, fluo, model, norm, iter_tr)
  res <- analyze_qpcR(all_fits, type, takeoff)
  res <- cbind(res, deltaF = calc_deltaF(input, cyc, fluo))
  res
})