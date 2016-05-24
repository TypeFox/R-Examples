smooth.fdPar <- function(fdobj, Lfdobj=NULL,
         lambda=0, estimate=TRUE, penmat=NULL){
##
## 1.  fdPar
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 2.  smooth.fd
##
  smooth.fd(fdobj, fdP)
}
