.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'InterVA4' package as:\n",
         "Zehang R. Li, Tyler H. McCormick and Samuel J. Clark (2014). ",
         "InterVA4: An R package to analyze verbal autopsy data. ",
         "Working paper no. 146, Center for Statistics and the Social Sciences, University of Washington\n"),
      domain = NULL,  appendLF = TRUE )
}