##' @export

###############################################################################
`gtExEs.default` <- function(proc){
ExePath <- system.file("RPRL", package = "MultiSV")
Exesrp <- (proc)
Exec <- file.path(ExePath, Exesrp)
Pexe <- Sys.which("perl")
paste(Pexe,Exec,sep=" ")
}


