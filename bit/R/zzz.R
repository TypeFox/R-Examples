.onLoad <- function(lib, pkg) {
  ##library.dynam("bit", pkg, lib) use useDynLib(bit) in NAMESPACE instead
  ##packageStartupMessage("Loading package bit ", packageDescription("bit", fields="Version"))
  bit_init()
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Attaching package bit")
  packageStartupMessage("package:bit (c) 2008-2012 Jens Oehlschlaegel (GPL-2)")
  packageStartupMessage("creators: bit bitwhich")
  packageStartupMessage("coercion: as.logical as.integer as.bit as.bitwhich which")
  packageStartupMessage("operator: ! & | xor != ==")
  packageStartupMessage("querying: print length any all min max range sum summary")
  packageStartupMessage("bit access: length<- [ [<- [[ [[<-")
  packageStartupMessage("for more help type ?bit")
}

.onUnload <- function(libpath){
   packageStartupMessage("Unloading package bit")
   bit_done()
   library.dynam.unload("bit", libpath)
}
