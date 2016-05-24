.onLoad <- function(lib, pkg) {
  if("catnet"%in%loadedNamespaces())
    unloadNamespace("catnet")
  dotpath <- as.character(Sys.getenv("R_DOTVIEWER"))
  ##if(dotpath == "")
  ##  packageStartupMessage("No plotting capabilities detected. Type 'help(sdnet)'\n")
}

