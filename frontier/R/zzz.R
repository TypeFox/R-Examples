.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'frontier' package as:\n",
         "Tim Coelli and Arne Henningsen (2013). ",
         "frontier: Stochastic Frontier Analysis. ",
         "R package version 1.0. ",
         "http://CRAN.R-Project.org/package=frontier.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'frontier' package, ",
         "please use a forum or 'tracker' at frontier's R-Forge site:\n",
         "https://r-forge.r-project.org/projects/frontier/"),
      domain = NULL,  appendLF = TRUE )
}
