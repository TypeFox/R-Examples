.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'InSilicoVA' package as:\n",
         "Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C. Crampin, Kathleen Kahn and Samuel J. Clark (2014). ",
         "Probabilistic cause-of-death assignment using verbal autopsies, ",
         "Journal of the American Statistical Association, to appear\n"),
      domain = NULL,  appendLF = TRUE )
}