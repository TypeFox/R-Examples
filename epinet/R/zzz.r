######################################################################
#
# zzz.r
#
# .First.lib is run when the package is loaded with library(epinet)
#
######################################################################

.onAttach <- function(lib, pkg){
  info <- packageDescription("epinet")
  packageStartupMessage(
    paste('\nepinet: version ', info$Version, ', created on ', info$Date, '\n', sep="")
    )
}
