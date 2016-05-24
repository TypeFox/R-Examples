#
# Do not use .Last.lib [CRAN check (3.0.0) causes note]
#
.onAttach <- function(...) { 
  pack <- "kitagawa"
  packageStartupMessage(
    sprintf("Loaded %s (%s) -- Spectral response of water wells", pack, utils::packageVersion(pack)))
}
