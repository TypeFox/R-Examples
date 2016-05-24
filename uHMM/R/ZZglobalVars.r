### Pour ouverture automatique de l'interface lors du chargement du package
#.onLoad <- function(libname = find.package("uHMM"), pkgname = "uHMM") {
#  uHMMinterface()
#}

# Pour eviter notes "<anonymous>: no visible binding for global variable"
globalVariables(c("CstateSeq","CsymbolSeq","selectedNames","estimatedHMM","symbCentersNorm","symbCenters","normParams","uHMMenv",
                  "nbSymbols","gap","MstateSeq","trainingRows","MarelCarnot","validationPeriod"
                  ))