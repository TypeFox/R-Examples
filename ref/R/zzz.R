.onAttach <- function(libname, pkgname){
  packageStartupMessage("dim(refdata) and dimnames(refdata) no longer allow parameter ref=TRUE, use dim(derefdata(refdata)), dimnames(derefdata(refdata)) instead")
}