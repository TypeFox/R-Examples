#' W' balance.
#'
#' Generate a vector of W' balance values from time and power data. The
#' underlying algorithm is published in Skiba \emph{et al.} (2012). \code{Wbal}
#' is a wrapper for the Rcpp function \code{Wbal_}.
#'
#' The algorithm used here, while based on Dr Phil Skiba's model, differs in
#' that values are positive as opposed to negative. The original published model
#' expressed W' balance as W' minus W' expended, the latter recovering with an
#' exponential time course when P < CP. An issue with this approach is that an
#' athlete might be seen to go into negative W' balance. Hence, to avoid
#' assumptions regarding available W', this algorithm returns W' expended (and
#' its recovery) as positive values; i.e. a ride is begun at 0 W' expended, and
#' it will \emph{increase} in response to supra-CP efforts.
#'
#' It is advisable on physiological grounds to enter smoothed power values to
#' the function, hence this is the default behaviour. If nothing else, this
#' prevents an unrealistic inflation of W' values that are inconsistent with
#' estimates derived from power-time modelling.
#'
#' The essence of the algorithm can be seen in the function
#' \href{https://github.com/jmackie4/cycleRtools/blob/master/tests/testthat/test_Wbal.R}{test
#' file.}
#'
#' Note that if there are \code{NA} values in the power column, these are
#' ignored and the correspoding W' expended value assumes that of the last
#' available power value. \code{NA} values are not allowed in the time column.
#'
#' @param data a data.frame/matrix object with time and power columns.
#' @param time character; name of the time (seconds) column in \code{data}.
#' @param pwr character; name of the power (watts) column in \code{data}.
#' @param t,P numeric vectors of time and power, respectively.
#' @param CP a critical power value for use in the calculation.
#' @param noisy logical; create smoother data by pooling power data into sub-
#'   and supra-CP sections.
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation.
#'
#' @examples \dontrun{
#' data(ridedata)
#'
#' ## Basic usage.
#' ridedata$Wexp.kJ <- Wbal(ridedata, timer.s, power.W, 310)
#'
#' ## Data can be noisy or "smooth"; e.g.
#' Wbal_noisy <- Wbal(ridedata, timer.s, power.W, 310, noisy = TRUE)
#' Wbal_smth  <- Wbal(ridedata, timer.s, power.W, 310, noisy = FALSE)
#'
#' ## Plot:
#' ylim <- rev(extendrange(Wbal_noisy))  # Reverse axes.
#'
#' plot(ridedata$timer.min, Wbal_noisy, type = "l", ylim = ylim,
#'      main = "NOISY")
#' plot(ridedata$timer.min, Wbal_smth, type = "l", ylim = ylim,
#'      main = "Smooooth")
#'
#' ## Example of NA handling.
#' d <- data.frame(t = seq_len(20), pwr = rnorm(20, 300, 50), Wexp.J = NA)
#' d[14:16, "pwr"] <- NA
#' d[, "Wexp.J"]   <- Wbal(d, "t", "pwr", CP = 290)
#' print(d)
#'
#' ## Using underlying Rcpp function:
#' Wbal_(t = 1:20, P = rnorm(20, 300, 50), CP = 300)  # Values are in joules.
#' }
#'
#' @return A numeric vector of W' balance values, in kilojoules or joules for
#'   \code{Wbal} or \code{Wbal_} respectively.
#'
#' @seealso \code{\link{plot.cycleRdata}}.
#'
#' @references Skiba, P. F., W. Chidnok, A. Vanhatalo, and A. M. Jones. Modeling
#'   the Expenditure and Reconstitution of Work Capacity above Critical Power.
#'   Med. Sci. Sports Exerc., Vol. 44, No. 8, pp. 1526-1532, 2012.
#'   \href{http://www.ncbi.nlm.nih.gov/pubmed/22382171}{PubMed link}.
#'
#' @export
Wbal <- function(data, time = "timer.s", pwr = "power.smooth.W",
                 CP = attr(data, "CP"), noisy = TRUE, character.only = FALSE) {

  if (!character.only) {
    time  <- data[, as.character(substitute(time))]
    pwr   <- data[, as.character(substitute(pwr))]
  }
  if (!noisy) {
    sec <- rep_len(0, length(pwr))
    sec[which(Diff(ifelse(pwr <= CP, 1, 0)) != 0) + 1] <- 1
    sec <- cumsum(sec)
    if (0 %in% sec) sec <- sec + 1
    pwr <- ave(pwr, sec, FUN = function(x) mean(x, na.rm = TRUE))
  }
  out <- Wbal_(time, pwr, CP)  # Rcpp.
  out <- out / 1000            # J to kJ.
  out
}
