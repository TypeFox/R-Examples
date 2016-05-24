
#   tnet R package :: http://toreopsahl.com/tnet/

.onAttach <- function(library, pkg) {
  packageStartupMessage("tnet: Analysis of Weighted, Two-mode, and Longitudinal networks.\nType ?tnet for help.")
}
