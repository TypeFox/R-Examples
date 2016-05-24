WeibullDiag <- function(formula, data = parent.frame(),
                        labels = names(m$strata))
{
  m <- survival::survfit(formula, data = data)
  if(!"strata" %in% names(m)) stop("No strata detected. No plot will be produced.")
  cp <- 1:length(m$strata)
  plot(m, fun = "cloglog",
       col = cp, lty = cp, lwd = 2, mark.time = FALSE,
       xlab = "Survival Time (log scale)", ylab = "Log Cumulative Hazard",
       main = "Weibull Regression\nDiagnostic Plot"
  )
  legend("topleft", legend = labels, col = cp,
         lwd = 2, lty = cp,
         bty = "n")
  pCH <- list(x = log(m$time), y = log(-log(m$surv)),
              strata = rep(names(m$strata), m$strata))
}


