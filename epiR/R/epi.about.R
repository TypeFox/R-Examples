"epi.about" <- function()
{
  cat("\n")
  cat("-----------------------------------------------------------\n")
  ver <- packageDescription("epiR", lib.loc = NULL, fields = "Version")
  cat(paste("epiR version", ver))
  cat("\n")
  cat("Tools for the Analysis of Epidemiological Data")
  cat("\n")
  cat("See http://fvas.unimelb.edu.au/veam for details.")
  cat("\n")
  cat("-----------------------------------------------------------\n")
  invisible()
}
