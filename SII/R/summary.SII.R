`summary.SII` <-
function (object, digits = 2, ...) 
{
  cat("\n")
  cat("Speech Intellibility Index\n")
  cat("--------------------------\n")
  cat("\n")
  cat("Call:\n")
  cat("\n")
  print(object$call)
  cat("\n")
  methodStr <- switch(object$method,
                      "interpolate"="Interpolation + Critical Band",
                      "critical"="Critical band SII procedure",
                      "equal-contributing"="Equally-contrlbutlng (17 band) critical-band SII procedure ",
                      "one-third octave"="One-third octave band SII procedure",
                      "octave"="Octave band SII procedure"
                      )
                      
  cat("Method:", methodStr, "\n" )
  cat("\n")
  cat("Calculation table:\n")
  cat("\n")
  print(object$table, digits=digits)
  cat("\n")
  cat("SII:", round(object$sii, digits), "\n")
  cat("\n")
  invisible(object)
}

