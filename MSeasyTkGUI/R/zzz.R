.onAttach <- function(lib,pkg)
{
  ###library(Rcmdr)
  #tkdestroy(commanderWindow)

  MSeasyTkGUI()
  packageStartupMessage("initializing ...", appendLF = FALSE)
}


