#' Calculates the thresholds
#'
#' \code{thresholds} alculates the thresholds
#' @param qp output from quickpsy
#' @param prob Probability to calculate the threshold.
#' @param log Use \code{TRUE}, if the logarithm of the independent variable
#' has been used to fit the curves (default is \code{FALSE}).
#' @export
#' @examples
#' library(MPDiR) # contains the Vernier data
#' data(Vernier) # ?Venier for the reference
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B =20,
#'                 thresholds = FALSE)
#' thresholds(fit, prob = .5)
#' @export
thresholds <- function(qp, prob = NULL, log = FALSE) {
  if (is.null(prob)) stop('You need to specify the value of prob', call. = FALSE)
    qp$par %>% do(one_threshold(., prob, log, qp$groups,
                               qp$funname, qp$guess, qp$lapses, qp$curves))
}



