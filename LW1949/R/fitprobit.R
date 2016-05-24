#' Fit a Probit Regression to Dose-Effect Data
#'
#' Fit a probit regression to dose-effect data, using the log10 of the dose as
#'   the response.
#' @param dat
#'   A data frame of toxicity data, including at least three variables:
#'     dose (the concentration of the tested chemical),
#'     ntot (the number of individuals tested), and
#'     nfx (the number of affected individuals).
#' @return
#'   A an object of class \code{\link{glm}}.
#' @export
#' @import
#'   stats
#' @details
#'   Only those rows with \code{dose > 0}, \code{ntot > 0}, and \code{nfx >= 0}
#'   are used in fitting the model.
#' @examples
#' toxdat <- data.frame(
#'  dose=c(0.05, 0.0625, 0.125, 0.25, 0.5, 1),
#'   ntot=rep(8, 6),
#'   nfx = c(0, 1, 4, 4, 6, 8))
#' fitprobit(toxdat)

fitprobit <- function(dat) {
  if (!is.data.frame(dat)) stop("Input must be a data frame.")
  if (any(is.na(match(c("dose", "ntot", "nfx"), names(dat))))) {
    stop("Input must include at least three variables: dose, ntot, nfx.")
  }
  sel <- with(dat, !is.na(dose) & dose>0 & !is.na(ntot) & ntot>0 &
      !is.na(nfx) & nfx>=0)
  if (sum(sel) < 1) stop("Data frame contains no rows of valid data.")
  glm(cbind(nfx, ntot-nfx) ~ log10(dose), family=binomial(link=probit),
    data=dat[sel, ])
  }
