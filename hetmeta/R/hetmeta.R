#' Deriving Measures Of Heterogeneity
#'
#' @description
#' The "\code{hetmeta}" implements the most common measures of heterogenity in meta-analysis.
#'
#' @param model an object of class "\code{\link{rma.uni}}".
#'
#' @details
#' The "\code{hetmeta}" function calculates estimates for several heterogeneity measures in
#' meta-analysis based on a meta-analytic model of class \code{\link{rma.uni}}
#' (see \code{\link{metafor-package}} for more details).
#'
#' Specifically, the measures derived in the function are the $R_b$, $I^2$, and $R_I$.
#' To complement those measures, the Dersimonian-Laird $Q$ test is presented, together with
#' the coefficient of variation of the pooled estimate $CV_b$, coefficient of variation
#' of the within-study variances, and the typical within-variance terms as
#' defined in the $I^2$ and $R_I$. See references for more details.
#'
#' @return
#' The \code{hetmeta} function returns an object of class "\code{hetmeta}" as described in \code{\link{hetmetaObject}}.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#'
#' @references
#'
#' Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. A new measure of between-studies
#' heterogeneity in meta-analysis. 2016. \emph{Stat. Med.} In Press.
#'
#' Takkouche B, Khudyakov P, Costa-Bouzas J, Spiegelman D. Confidence Intervals for Heterogeneity Measures in Meta-analysis. \emph{Am. J. Epidemiol.} 2013:kwt060.
#'
#' Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Stat. Med.} 2002; 21(11):1539-1558.
#'
#' Takkouche B, Cadarso-Suarez C, Spiegelman D. Evaluation of old and new tests of heterogeneity in epidemiologic meta-analysis. Am. J. Epi- demiol. 1999; 150(2):206-215.
#'
#' @seealso \code{\link{hetmeta-package}}, \code{\link{metafor}}
#'
#' @export hetmeta
#'
#' @examples
#' ## load data
#' dat <- get(data(dat.gibson2002))
#'
#' ## random-effects model analysis of the standardized mean differences
#' dat <- escalc(measure = "SMD", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i,
#'               sd2i = sd2i, n2i = n2i, data = dat)
#' res <- rma(yi, vi, data = dat, method = "REML")
#'
#' ## heterogeneity measures
#' hetmeta(res)
#'
#'
#' ## load BCG vaccine data
#' data(dat.bcg)
#'
#' ## random-effects model of log relative risks
#' dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#' res <- rma(yi, vi, data=dat)
#'
#' ## heterogeneity measures
#' hetmeta(res)
#'

hetmeta <- function(model){

   if (!is.element("rma.uni", class(model)))
      stop("Argument 'x' must be an object of class \"rma.uni\".")
   k <- model$k
   tau2 <- model$tau2
   wfi <- 1/model$vi
   si <- function(n) sum(wfi^n)
   Q <- model$QE

   # Typical within-study variance (for I^2 and R_I)
   s2_I2 <- (k-1)*si(1)/(si(1)^2 - si(2))
   s2_Ri <- k/si(1)
   # Estimate of the (squared of) the coefficient of variation of (fixed-effects) weights.
   cv2 <- k*si(2)/(si(1)^2) - 1
   # Estimate of the coefficient of variation of of within stydy variances.
   cv_vi <- k*si(2)/(si(1)^2) - 1
   # Measures of heterogeneity
   #I2 <- model$I2
   Ri <- 100*tau2/(tau2 + k/si(1))
   Rb <- 100*tau2/(k*model$se^2)
   CVb <- sqrt(tau2/c(model$b^2))

   # Approximated asymptotic stderr
   se_Q <- sqrt(2*(k - 1) + 4*(si(1) - si(2)/si(1))*tau2 +
                   2*(si(2) - 2*si(3)/si(1) + si(2)^2/si(1)^2)*tau2^2)
   se_I2 <- sqrt(se_Q^2 * (k - 1)^2/Q^4)
   se_Ri <- sqrt(se_Q^2 * (k - 1 - cv2)^2/(Q - cv2)^4)
   # 'as' need for the variance of R_I
   as <- model$vi * (si(1) - si(2)/si(1))
   se_Rb <- sqrt(se_Q^2 * (sum(as/(Q + as - (k-1))^2)/k)^2)
   se_CVb <- sqrt(tau2 * c(model$se^2)/c(model$b^4) + model$se.tau2^2/
                               (4*c(model$b^2)*tau2))

   # Results
   model <- c(model, list(Rb = Rb, Ri = Ri, CVb = CVb,
                          se_I2 = se_I2, se_Ri = se_Ri, se_Rb = se_Rb, se_CVb = se_CVb,
                          s2_I2 = s2_I2, s2_Ri = s2_Ri, cv_vi = cv_vi))

   class(model) <- c("hetmeta", "rma.uni")
   model
}
