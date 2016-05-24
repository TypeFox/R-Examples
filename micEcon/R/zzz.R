.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nIf you have questions, suggestions, or comments ",
         "regarding one of the 'micEcon' packages, ",
         "please use a forum or 'tracker' at micEcon's R-Forge site:\n",
         "https://r-forge.r-project.org/projects/micecon/"),
      domain = NULL,  appendLF = TRUE )
}
