.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE) 
#    require(setRNG)
}

setClassUnion("vectororNULL", c("vector","NULL"))

.distrSimoptions <- list(
                      MaxNumberofPlottedObs = 4000,
                      MaxNumberofPlottedObsDims = 6,
                      MaxNumberofPlottedRuns = 6,
                      MaxNumberofSummarizedObsDims = 6,
                      MaxNumberofSummarizedRuns = 6,
                      MaxNumberofPrintedObsDims = 10,
                      MaxNumberofPrintedObs = 10,
                      MaxNumberofPrintedRuns = 100
                      )
  

.onAttach <- function(library, pkg)
{
  unlockBinding(".distrSimoptions", asNamespace("distrSim"))
    msga <- gettext(
    "Some functions from package 'stats' are intentionally masked ---see distrSimMASK().\n"
                   )
  buildStartupMessage(pkg="distrSim", msga, packageHelp=TRUE, library=library, 
               #     MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
  VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\")."))
###
  invisible()
}

distrSimMASK <- function(library = NULL) 
{
    infoShow(pkg = "distrSim", filename="MASKING", library = library)
}
