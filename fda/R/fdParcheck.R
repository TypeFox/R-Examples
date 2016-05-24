fdParcheck = function (fdParobj) {
  if (!inherits(fdParobj, "fdPar")) {
    if (inherits(fdParobj, "fd") || inherits(fdParobj, "basisfd")) {
        fdParobj <- fdPar(fdParobj)
    } else
        stop(paste("'fdParobj' is not a functional parameter object,",
               "not a functional data object, and",
               "not a basis object."))
  }

return(fdParobj)

}
