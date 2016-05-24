.onLoad<-function(lib, pkg) {
    library.dynam("isocir", pkg, lib)
    ## echo output to screen
    #packageStartupMessage("initializing ...", appendLF = FALSE)
    #packageStartupMessage(" ISOtonic Inference for CIRcular Data\n")
    ## assuming you would need the packages combinat and circular.
    #require(combinat)
    #require(circular)
    #packageStartupMessage(" done")
}