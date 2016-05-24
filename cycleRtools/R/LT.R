#' Lactate Thresholds
#'
#' Model lactate threshold markers from work rate (power) and blood lactate
#' values. Requires package "pspline".
#'
#' @param WR a numeric vector of work rate values. Typically these would be the
#'   work rates associated with stages in an incremental exercise test.
#' @param La a numeric vector of blood lactate values (mmol/L) associated with
#'   the stages described in \code{WR}.
#' @param sig_rise numeric; a rise in blood [Lactate] that is deemed
#'   significant. Default is 1.5 mmol/L.
#' @param plots should outputs be plotted?
#'
#' @details This function is a slightly modified version of that written by
#'   Newell \emph{et al.} (2007) and published in the Journal of Sport Sciences
#'   (see references). The original source code, which also includes other
#'   functions for lactate analysis, can be found
#'   \href{http://www.nuigalway.ie/maths/jn/Lactate/html/lactate-r.html}{here}.
#'
#' @references John Newell , David Higgins , Niall Madden , James Cruickshank ,
#'   Jochen Einbeck , Kenny McMillan & Roddy McDonald (2007) Software for
#'   calculating blood lactate endurance markers, Journal of Sports Sciences,
#'   25:12, 1403-1409, \href{http://dx.doi.org/10.1080/02640410601128922}{DOI}.
#'
#' @return a data frame of model outputs, and optionally a matrix of plots.
#'
#' @examples
#' # This data is included with Newell et al's source code.
#' WR <- c(50, 75, 100, 125, 150, 175, 200, 225, 250)
#' La <- c(2.8, 2.4, 2.4, 2.9, 3.1, 4.0, 5.8, 9.3, 12.2)
#' LT(WR, La, 1.5, TRUE)
#'
#' @seealso Newell \emph{et al.}'s
#'   \href{https://orreco.shinyapps.io/lactate/}{Shiny app.}
#'
#' @export
LT <- function(WR, La, sig_rise = 1.5, plots = TRUE) {
  if (length(WR) != length(La) | any(is.na(c(WR, La))))
    stop("Invalid inputs", call. = FALSE)
  # Functions
  d2lmax.discrete <- function(WR, La) {
    mt2 <- numeric()
    for (t in 1:length(La))
      mt2[t] <- La[t + 2] - 2 * La[t + 1] + La[t]
    d2lmax.disc <- WR[mt2 == max(mt2, na.rm = T)][1]
    ifelse(length(WR) <= 4, NA, d2lmax.disc)
  }
  peaks <- function(series, span = 3) {
    z <- embed(series, span)
    max.col(z) == 1 + span %/% 2
  }
  LT.est <- function(La2, WR2) {
    LT <- round(optim(
      WR2[2], LT.bs.reg2, method = "BFGS", y = La2, x = WR2
    )$par, digits = 2)
    lWR <- log(WR2)
    lLa <- log(La2)
    LT2 <- round(optim(
      lWR[2], LT.bs.reg2, method = "BFGS", y = lLa, x = lWR
    )$par, 2)
    c(LT, exp(LT2))
  }
  LT.bs.reg2 <- function(break.pt, y, x) {
    lt.fit <-
      lm(y ~ lhs(x, break.pt) + rhs(x, break.pt), singular.ok = T)
    lt.mse <- (sum(lt.fit$residuals ^ 2)) / lt.fit$df.residual
    lt.mse
  }
  lhs <- function(x, break.pt) ifelse(x < break.pt, break.pt - x, 0)
  rhs <- function(x, break.pt) ifelse(x < break.pt, 0, x - break.pt)
  # D2Lmax bit------------------------------------------------------------------
  fit <-
    pspline::smooth.Pspline(WR, La, norder = 3, df = length(La) - 3, method = 2)
  La.grid <-
    seq(from = min(WR), to = max(WR), length = 1000)
  fit.deriv.2 <- predict(fit, La.grid, 2)
  cbind(La.grid, fit.deriv.2)
  d2lmax <- max(La.grid[peaks(fit.deriv.2, span = 3)])
  # D2Lmax Discrete-------------------------------------------------------------
  d2lmax2 <- d2lmax.discrete(WR, La)
  # DMax------------------------------------------------------------------------
  poly.coefs <- lm(La ~ WR + I(WR ^ 2) + I(WR ^ 3))$coefficients
  linbeta <-
    (La[length(La)] - La[1]) / (WR[length(WR)] - WR[1])
  # Need to solve where the first derivative of the poly.fit equals the slope of
  # the line joining the first and last La readings.
  dmax <-
    Re(polyroot(c(
      poly.coefs[2] - linbeta,2 * poly.coefs[3],3 * poly.coefs[4]
    )))
  dmax <- max(dmax[dmax <= max(WR)])
  # 4mmol estimate--------------------------------------------------------------
  est.4mmol <-
    Re(polyroot(c(
      poly.coefs[1] - 4,poly.coefs[2],poly.coefs[3],poly.coefs[4]
    )))
  est.4mmol <- max(est.4mmol[est.4mmol <= max(WR)])
  # 1mmom > baseline------------------------------------------------------------
  risepb <- La[1] + sig_rise
  Rise.mmol <-
    Re(polyroot(
      c(poly.coefs[1] - risepb, poly.coefs[2], poly.coefs[3], poly.coefs[4])
    ))
  Rise.mmol <- max(Rise.mmol[Rise.mmol <= max(WR)])
  # LT and Log log LT estimates-------------------------------------------------
  est.LT <- LT.est(La, WR)
  lt.fit <-
    lm(La ~ lhs(WR, est.LT[1]) + rhs(WR, est.LT[1]), singular.ok = T)
  lt.pred.y <-
    lt.fit$coef[1] + lt.fit$coef[2] * lhs(La.grid, est.LT[1]) +
    lt.fit$coef[3] * rhs(La.grid, est.LT[1])
  loglt.pred.y <-
    lt.fit$coef[1] + lt.fit$coef[2] * lhs(La.grid, est.LT[2]) +
    lt.fit$coef[3] * rhs(La.grid, est.LT[2])
  y.LT <-
    lt.fit$coef[1] + lt.fit$coef[2] * lhs(est.LT[1], est.LT[1]) +
    lt.fit$coef[3] * rhs(est.LT[1], est.LT[1])
  y.logLT <-
    lt.fit$coef[1] + lt.fit$coef[2] * lhs(est.LT[2], est.LT[2]) +
    lt.fit$coef[3] * rhs(est.LT[2], est.LT[2])
  # Some fits needed for plots.
  linear.fit <-
    lm(La[c(1, length(La))] ~ WR[c(1, length(La))])
  poly.fit <- lm(La ~ WR + I(WR ^ 2) + I(WR ^ 3))
  poly.predict <- predict(poly.fit, WR = La.grid)
  # Return----------------------------------------------------------------------
  Workrate <-
    round(c(est.LT, est.4mmol, Rise.mmol, dmax, d2lmax , d2lmax2), 2)
  out <- data.frame(
    Workrate,
    row.names = c(
      "LT", "LT.log.log", "4mmol", paste("Rise", sig_rise, "mmol"),
      "Dmax", "D2LMax", "D2LMax Discrete"
    )
  )
  if (!plots) return(out)
  # Plots-----------------------------------------------------------------------
  split.screen(c(2, 3))
  screen(1) #-------------------------------------------------------------------
  plot(WR, La, xlab = "Workrate", ylab = "La mmol/l", ylim = c(0, max(La) + 1),
       las = 1)
  title(main = paste("LT =", round(est.LT[1], 2)))
  abline(v = est.LT[1], lty = 2)
  screen(2) #-------------------------------------------------------------------
  plot(WR, La, xlab = "Workrate", ylab = "La mmol/l", ylim = c(0, max(La) + 1),
       las = 1)
  title(main = paste("LTloglog =", round(est.LT[2], 2)))
  abline(v = est.LT[2], lty = 2)
  screen(3) #-------------------------------------------------------------------
  plot(WR, La, xlab = "Workrate   ", ylab = "La mmol/l", ylim = c(0, max(La) + 1),
       las = 1)
  title(main = paste("Rise", sig_rise, "mmol =", round(Rise.mmol, 2)))
  abline(h = risepb, lty = 2)
  abline(v = Rise.mmol, lty = 2)
  screen(4) #-------------------------------------------------------------------
  plot(WR, La, xlab = "Workrate", ylab = "La mmol/l", ylim = c(0, max(La) + 1),
       las = 1)
  title(main = paste("FBLC 4 mmol =", round(est.4mmol, 2)))
  abline(h = 4, lty = 2)
  abline(v = est.4mmol, lty = 2)
  screen(5) #-------------------------------------------------------------------
  plot(WR, La, xlab = "Workrate", ylab = "La mmol/l", ylim = c(0, max(La) + 1),
       las = 1)
  title(main = paste("DMax =", round(dmax, 2)))
  abline(v = dmax, lty = 2)
  screen(6) #-------------------------------------------------------------------
  plot(La.grid, fit.deriv.2, xlab = "Workrate", ylab = "D2L/DW")
  abline(v = c(d2lmax), lty = 2)
  title(main = paste("D2LMax =", round(d2lmax, 2)))
  close.screen(all.screens = TRUE)
  #-----------------------------------------------------------------------------
  return(out)
}
