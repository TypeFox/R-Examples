.onLoad <- function(lib,pkg) {
    library.dynam("VBLPCM", pkg, lib)
    }

.onAttach <- function(lib,pkg) {
    packageStartupMessage("\n\n")
    packageStartupMessage(paste("\n", "Variational Bayes Latent Position Cluster Model for networks.","\n"),
        paste("VBLPCM", "Version", "2.4.3", 
        "Created on", "2014-01-29"), paste("\n Created and maintained by ", "Michael Salter-Townshend"), "\n")
    packageStartupMessage('For citation information type \'citation("VBLPCM")\'\n')
    packageStartupMessage('Type \'help(VBLPCM)\' to get started.\n')
    packageStartupMessage('Some worked examples are given by \'example(VBLPCM)\' \n')
} 

.onUnload <- function(libpath){
  library.dynam.unload("VBLPCM", libpath)
} 
