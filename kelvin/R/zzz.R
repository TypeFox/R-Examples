#
.onAttach <- function(...) { 
  pack <- 'kelvin'
  packageStartupMessage(
    sprintf("Loaded %s (%s) -- Solutions to the Kelvin differential equation", pack, utils::packageVersion(pack))
  )
}