#' Peak Test
#' 
#' Detect, separate and count positive and negative peaks, as well as peak-like
#' noise. Additionally, function calculates area of the peaks.
#' 
#' The localization of peaks is determined by the
#' \code{\link[pracma]{findpeaks}} function.  The area under the peak is
#' calculated by integration of approximating spline.
#' 
#' @name test_peaks
#' @aliases test_peaks test_peaks-methods test_peaks,adpcr-method
#' test_peaks,numeric-method test_peaks
#' @docType methods
#' @param x a vector containing the abscissa values (e.g., time, position) OR
#' an object of class \code{\linkS4class{adpcr}}.
#' @param y a vector of fluorescence value.
#' @param threshold a value, which defines the peak heights not to consider as
#' peak.
#' @param noise_cut a numeric value between 0 and 1. All data between 0 and
#' \code{noise_cut} quantile would be considered noise in the further analysis.
#' @param savgol logical value. If \code{TRUE}, Savitzky-Golay smoothing filter
#' is used.
#' @param norm logical value. If \code{TRUE}, data is normalised.
#' @param filter.q a vector of two numeric values. The first element represents
#' the quantile of the noise and the second one is the quantile of the negative
#' peaks.
#' @return A list of length 2. The first element is a data frame containing:
#' peak number, peak group (noise, negative, positive), position of the peak
#' maximum, area under the peak, peak width, peak height, position of the peak
#' and time resolution.
#' 
#' The second element contains smoothed data.
#' @author Stefan Roediger, Michal Burdukiewicz.
#' @references Savitzky, A., Golay, M.J.E., 1964. Smoothing and Differentiation
#' of Data by Simplified Least Squares Procedures. Anal. Chem. 36, 1627-1639.
#' @keywords smooth peak AUC noise
#' @export
#' @examples
#' 
#' 
#' data(many_peaks)
#' par(mfrow = c(3,1))
#' plot(many_peaks, type = "l", main = "Noisy raw data")
#' abline(h = 0.01, col = "red")
#' 
#' tmp.out <- test_peaks(many_peaks[, 1], many_peaks[, 2], threshold = 0.01, noise_cut = 0.1, 
#'                     savgol = TRUE)
#' plot(tmp.out[["data"]], type = "l", main = "Only smoothed")
#' abline(h = 0.01, col = "red")
#' abline(v = many_peaks[tmp.out[["peaks"]][, 3], 1], lty = "dashed")
#' 
#' tmp.out <- test_peaks(many_peaks[, 1], many_peaks[, 2], threshold = 0.01, noise_cut = 0.1, 
#'                     savgol = TRUE, norm = TRUE)
#' plot(tmp.out[["data"]], type = "l", main = "Smoothed and peaks detected")
#' abline(v = many_peaks[tmp.out[["peaks"]][, 3], 1], lty = "dashed")
#' for(i in 1:nrow(tmp.out$peaks)) {
#'   if(tmp.out$peaks[i, 2] == 1) {col = 1}
#'   if(tmp.out$peaks[i, 2] == 2) {col = 2}
#'   if(tmp.out$peaks[i, 2] == 3) {col = 3}
#'   points(tmp.out$peaks[i, 7], tmp.out$peaks[i, 6], col = col, pch = 19)
#' }
#' 
#' positive <- sum(tmp.out$peaks[, 2] == 3)
#' negative <- sum(tmp.out$peaks[, 2] == 2)
#' total <- positive + negative
#' 
#' 
NULL


test_peaks <- function (x, ...) {
  stop("Wrong class of 'x'", call. = TRUE, domain = NA)
}

setGeneric("test_peaks")


setMethod("test_peaks", signature(x = "numeric"), function(x, y, threshold = 0.05, 
                                                           noise_cut = 0.05, savgol = TRUE, 
                                                           norm = FALSE,
                                                           filter.q = c(0.7, 0.8)) {
  # Initial checking of input values
  if (is.null(x)) 
    stop("Enter 'x' value.", call. = TRUE, domain = NA)
  if (is.null(y)) 
    stop("Enter 'y' value.", call. = TRUE, domain = NA)
  if (length(x) != length(y)) 
    stop("'x' and 'y' differ in length.", call. = TRUE, domain = NA)
  AUCtest(x = x, y = y, threshold = threshold, noise_cut = noise_cut, savgol = savgol, 
          norm = norm, filter.q = filter.q)
  
})

setMethod("test_peaks", signature(x = "adpcr"), function(x, threshold = 0.05, 
                                                         noise_cut = 0.05, savgol = TRUE, 
                                                         norm = FALSE,
                                                         filter.q = c(0.7, 0.8)) {
  # Initial checking of input values
  if (slot(x, "type") != "fluo") 
    stop("'adpcr' object must have type 'fluo'.", call. = TRUE, domain = NA)
  x_vals <- slot(x, ".Data")[[1]]
  y_vals <- slot(x, ".Data")[[2]]
  AUCtest(x = x_vals, y = y_vals, threshold = threshold, noise_cut = noise_cut, savgol = savgol, 
          norm = norm, filter.q = filter.q)
  
})