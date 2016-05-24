
.onAttach <- function(libname,pkgname){
#    msg <- ""
#    nfc <- "Use suppressPackageStartupMessages(library(editrules)) to suppress this message on loading editrules."
#    packageStartupMessage(msg)
#    packageStartupMessage(nfc)

# default lpSolveAPI options
options(er.lpcontrol =   
  list( 
    presolve = "rows"    # move univariate constraints into bounds
    , epsint = 1e-15
#   , epssel = 1e-15
#   , epsb = 1e-15
#   , epsd = 1e-15
    , epspivot = 1e-15
   )
)
}


