.onLoad <- function (lib, pkg) {
  library.dynam("Rrdrand", pkg, lib)
  if(hasRDRAND()==TRUE){
    RNGkind("user-supplied")
  }

}
.onAttach <-function (lib, pkg) {
  if(hasRDRAND()==FALSE){
    packageStartupMessage("WARNING : This CPU does not support RDRAND:  RNGkind has not been changed.")
  }
}

.onUnload <- function (libpath) {
  RNGkind("default")
  library.dynam.unload("Rrdrand", libpath)
}
