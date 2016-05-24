print.ecm <- function(x, ...)
{
  if (inherits(x, "ecmSymFit")) {
      ti = "ECM - Symmetric + linear cointegration - "
  } else {
      if (inherits(x, "ecmAsyFit") & x$model == "linear") {
        ti = "ECM - Asymmetric + linear cointegration - "
      } else {
        ti = "ECM - Asymmetric + nonlinear threshold cointegration - "
      }
  }

  cat("\n===============================================================")
  cat(paste("\n", ti, " \"", x$name.x, "\"", "\n", sep = ""))       
  cat("===============================================================\n")       
  print(summary(x$ecm.x))
  
  cat("\n===============================================================")
  cat(paste("\n", ti, " \"", x$name.y, "\"", "\n", sep = ""))       
  cat("===============================================================\n")       
  print(summary(x$ecm.y))
} 