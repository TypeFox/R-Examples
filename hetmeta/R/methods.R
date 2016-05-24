#' hetmeta Object
#'
#' @description An object returned by \code{hetmeta} function, inheriting from class "\code{rma.uni}"
#'
#' @name hetmetaObject
#' @docType class
#'
#' @details
#'
#' An object of class "\code{hetmeta}". The object is derived from an object of class \code{\link{rma.uni}}.
#' In addition to thatm it has the following components:
#' Objects of class "\code{hetmeta}" are lists with defined components.
#' \tabular{ll}{
#' \code{Rb} \tab value of \eqn{R_b}{Rb}, which quantifies the proportion of the between-study heterogeneity
#' relative to the variance of the pooled random effects estimate. \cr
#' \code{Ri} \tab value of \eqn{R_I}, whihc quantifies the proportion of the variance of the effect
#'  estimate due to between-studies variation. \cr
#' \code{CVb} \tab value of \eqn{CV_b}, the between-studies coefficient of variation.  \cr
#' \code{se_Rb} \tab the sandard error of \eqn{R_b2} derived using the delta method. \cr
#' \code{se_I2} \tab the sandard error of \eqn{I^2} derived using the delta method. \cr
#' \code{se_Ri} \tab the sandard error of \eqn{R_I} derived using the delta method. \cr
#' \code{se_CVb} \tab the sandard error of \eqn{CV_b} derived using the delta method. \cr
#' \code{s2_I2} \tab the "typical" within-study variance as defined in the \eqn{I^2} \cr
#' \code{s2_Ri} \tab the "typical" within-study variance as defined in the \eqn{R_I} \cr
#' \code{cv_vi} \tab value of the coefficient of variation of the within-study variances.  \cr
#' }
#'
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#'
#' @seealso \code{\link{hetmeta}}, \code{\link{hetmeta-package}}

NULL

#' Printing hetmeta Results
#'
#' @description Print function for objects of class "\code{hetmeta}".
#'
#' @param x an object of class \code{hetmeta} produced by \code{\link{hetmeta}}.
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#'
#' @seealso \code{\link{hetmeta}}
#'
#' @rdname print.hetmeta
#' @method print hetmeta
#' @export
#'
#' @examples
#' #To be included

print.hetmeta <- function(x, digits, ...){
   if (!is.element("hetmeta", class(x)))
      stop("Argument 'x' must be an object of class \"hetmeta\".")
   if (missing(digits))
      digits <- x$digits
   cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
   ncutoff <- as.numeric(cutoff)

   if (x$method == "FE"){
      if (x$int.only) {
         cat("Fixed-Effects Model (k = ", x$k, ")", sep = "")
      } else {
         cat("Fixed-Effects with Moderators Model (k = ",
             x$k, ")", sep = "")
      }
   } else {
      if (x$int.only) {
         cat("Random-Effects Model (k = ", x$k, "; ", sep = "")
      }
      else {
         cat("Mixed-Effects Model (k = ", x$k, "; ", sep = "")
      }
      cat("tau^2 estimator: ", x$method, ")", sep = "")
   }
   cat("\n\n")

   if (x$method != "FE") {
      if (x$int.only) {
         cat("tau^2 (estimated amount of total heterogeneity):  ",
             formatC(x$tau2, digits = ifelse(abs(x$tau2) <=
                                                .Machine$double.eps * 10, 0, digits), format = "f"),
             ifelse(is.na(x$se.tau2), "", paste0(" (SE = ",
                                                 formatC(x$se.tau2, digits = digits, format = "f"),
                                                 ")")), "\n", sep = "")
      } else {
         cat("tau^2 (estimated amount of residual heterogeneity):  ",
             formatC(x$tau2, digits = ifelse(abs(x$tau2) <=
                                                .Machine$double.eps * 10, 0, digits), format = "f"),
             ifelse(is.na(x$se.tau2), "", paste0(" (SE = ",
                                                 formatC(x$se.tau2, digits = digits, format = "f"),
                                                 ")")), "\n", sep = "")
      }
      cat("tau (square root of estimated tau^2 value):       ",
          ifelse(x$tau2 >= 0, formatC(sqrt(x$tau2), digits = ifelse(x$tau2 <=
                                                                       .Machine$double.eps * 10, 0, digits), format = "f"),
                 NA), "\n", sep = "")

      ## Measures of heterogeneity
      cat("H^2 (total variability / sampling variability):   ",
          ifelse(is.na(x$H2), NA, formatC(x$H2, digits = max(0, digits),
                                          format = "f")), "\n", "\n", sep = "")
      if (!is.null(x$R2))
         cat("R^2 (amount of heterogeneity accounted for):            ",
             ifelse(is.na(x$R2), NA, formatC(x$R2, digits = max(0, digits),
                                             format = "f")), "%", "\n", sep = "")

      cat("R_b:     ",
          formatC(x$Rb, digits = max(0, digits-2), format = "f"), "%",
          paste0(" (SE = ",
                 formatC(100*x$se_Rb, digits = digits, format = "f"),
                 ")"), "\n", sep = "")
      cat("I^2:     ",
          ifelse(is.na(x$I2), NA, formatC(x$I2, digits = max(0, digits-2),
                                          format = "f")), "%",
          paste0(" (SE = ",
                 formatC(100*x$se_I2, digits = digits, format = "f"),
                 ")"), "\n", sep = "")
      cat("R_I:     ",
          formatC(x$Ri, digits = max(0, digits-2), format = "f"), "%",
          paste0(" (SE = ",
                 formatC(100*x$se_Ri, digits = digits, format = "f"),
                 ")"), "\n", sep = "")
      cat("\n")
      cat("CVb between-study coefficient of variation:             ",
          formatC(x$CVb, digits = digits, format = "f"),
          paste0(" (SE = ",
                 formatC(x$se_CVb, digits = digits, format = "f"),
                 ")"), "\n", sep = "")
   }

   cat("coefficient of variation of the within-study variances: ",
       formatC(sqrt(x$cv_vi), digits = digits, format = "f"), "\n", sep = "")
   cat("'typical' within-study variance defined in I^2:         ",
       formatC(x$s2_I2, digits = digits, format = "f"), "\n", sep = "")
   cat("'typical' within-study variance defined in R_I:         ",
       formatC(x$s2_Ri, digits = digits, format = "f"), sep = "")
   cat("\n\n")

   if (!is.na(x$QE)) {
      QEp <- x$QEp
      if (QEp > ncutoff) {
         QEp <- paste("=", formatC(QEp, digits = digits, format = "f"))
      }
      else {
         QEp <- paste0("< ", cutoff)
      }
      if (x$int.only) {
         cat("Test for Heterogeneity: \n")
         cat("Q(df = ", x$k - x$p, ") = ", formatC(x$QE, digits = digits,
                                                   format = "f"), ", p-val ", QEp, "\n\n", sep = "")
      }
      else {
         cat("Test for Residual Heterogeneity: \n")
         cat("QE(df = ", x$k - x$p, ") = ", formatC(x$QE,
                                                    digits = digits, format = "f"), ", p-val ", QEp,
             "\n\n", sep = "")
      }
   }
   QMp <- x$QMp
   if (QMp > ncutoff) {
      QMp <- paste("=", formatC(QMp, digits = digits, format = "f"))
   }
   else {
      QMp <- paste0("< ", cutoff)
   }

   invisible()
}


#' Confidence Intervals for 'hetmeta' Objects
#'
#' @description The function calculates confidence intervals for the heterogeneity measures in a 'hetmeta' object.
#'
#' @param object an object of class \code{hetmeta} produced by \code{\link{hetmeta}}.
#' @param parm this argument is here for compatability with the generic function confint, but is (currently) ignored.
#' @param level numerical value between 0 and 100 specifying the confidence interval level (if unspecified, the default is to take the value from the object).
#' @param digits an integer specifying the number of digits to which printed results must be rounded.
#' @param \dots further arguments passed to or from other methods.
#'
#' @details The confidence intervals are constructed based on the (asymptotic) normal
#' distribution of the estimators. Standard error are derived using the delta method.
#' See the references for more details.
#'
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#'
#' @seealso \code{\link{hetmeta}}
#'
#' @references
#'
#' Takkouche B, Khudyakov P, Costa-Bouzas J, Spiegelman D. Confidence Intervals for Heterogeneity Measures in Meta-analysis. \emph{Am. J. Epidemiol.} 2013:kwt060.
#'
#' Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. A new measure of between-studies heterogeneity in meta-analysis. 2016. \emph{Stat. Med.} In Press.
#'
#' @rdname confint.hetmeta
#' @method confint hetmeta
#' @export
#'
#' @examples
#' ## load BCG vaccine data
#' data(dat.bcg)
#'
#' ## random-effects model of log relative risks
#' dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#' res <- rma(yi, vi, data=dat)
#'
#' ## heterogeneity measures
#' het <- hetmeta(res)
#' confint(het)

confint.hetmeta <- function(object, parm, level, digits, ...){
   x <- object
   if (!is.element("hetmeta", class(x)))
      stop("Argument 'object' must be an object of class \"hetmeta\".")
   if (x$method == "FE")
      stop("Method' available only for random-effect models.")
   if (missing(level))
      level <- x$level
   if (missing(digits))
      digits <- x$digits
   alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)

   tab <- cbind(stat <- c(x$Rb, x$I2, x$Ri, x$CVb),
                se.stat <- c(100*x$se_Rb, 100*x$se_I2, 100*x$se_Ri, x$se_CVb),
                ci.lb = stat - qnorm(1-alpha/2)*se.stat,
                ci.ub = stat + qnorm(1-alpha/2)*se.stat)
   tab[-4, 3] <- pmax(0, tab[-4, 3])
   tab[-4, 4] <- pmin(100, tab[-4, 4])
   colnames(tab) <- c("estimate", "se", "ci.lb", "ci.ub")
   rownames(tab) <- c("R_b (%)" ,"I^2 (%)", "R_I (%)", "CV_b")
   table <- formatC(tab, digits = digits, format = "f")

   print(table, quote = FALSE, right = TRUE, print.gap = 2)
   cat("\n")
   invisible()
}
