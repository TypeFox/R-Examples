.onAttach <- function(...) {
  if (stats::runif(1) > 0.15) return()
  else (interactive())
  {
    packageStartupMessage("Type 'help(package=dkDNA)' for the package overview")
  }
}
