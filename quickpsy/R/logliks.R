#' Calculates the loglikelihoods
#'
#' \code{logliks} calculates the loglikelihoods.
#' @param qp output from quickpsy
#' @export
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 20)
#' logliks(fit)
#' @export
logliks <- function(qp) {
  qp$averages %>%
    do(one_loglik(., qp$x, qp$k, qp$n, qp$psyfunguesslapses, qp$groups,
                  qp$par))
}


