#' Predicted probabilities
#'
#' \code{ypred} calculates the predicted probabilities at the values of the
#' explanatory variable.
#' @param qp output from quickpsy
#' @examples
#' library(MPDiR) # contains the Vernier data
#' data(Vernier) # ?Venier for the reference
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 20)
#' ypred(fit)
#' @export
ypred <- function(qp) {
  qp$par %>% do(one_ypred(., log, qp$groups, qp$averages, qp$x,
                           qp$psyfunguesslapses))

}
