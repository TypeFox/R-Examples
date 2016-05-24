.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'maxLik' package as:\n",
         "Henningsen, Arne and Toomet, Ott (2011). ",
         "maxLik: A package for maximum likelihood estimation in R. ",
         "Computational Statistics 26(3), 443-458. ",
         "DOI 10.1007/s00180-010-0217-1.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'maxLik' package, ",
         "please use a forum or 'tracker' at maxLik's R-Forge site:\n",
         "https://r-forge.r-project.org/projects/maxlik/"),
      domain = NULL,  appendLF = TRUE )
}
