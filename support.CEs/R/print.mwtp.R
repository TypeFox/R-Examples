print.mwtp <- function(x, digits = max(3, getOption("digits") - 3), scientific = FALSE, ...)
{
# Name: print.mwtp
# Title: print() for S3 class "mwtp"
# Arguments:
#  x            An object of S3 class "gofm".
#  digits       The number of significant digits.
#  scientific   Scores are encoded in scientific format.
#  ...          Arguments passed to format().



# display table of MWTPs

  cat("\n")
  print(format(x$mwtp.table, digits = digits, scientific = scientific, ...), 
        quote = FALSE, right = TRUE)
  cat("\nmethod =", ifelse(x$method == "kr", "Krinsky and Robb", "Delta"), "\n")
  cat("\n")

  invisible(x)
}

