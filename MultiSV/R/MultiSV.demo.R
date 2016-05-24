##' Function with MultiSV input example
##'
##' 
##' 
##'
##' 
##' 
##' @example R/MultiSV.example.R
##' @author Khurram \enc{Maqbool}{Maqbool}
##' @export
GetTestCfgFile <- function () {
ExtDtPth <- system.file("extdata", package = "MultiSV")
CfgFileName <- ("MultiSV.conf")
CfgFile <- file.path(ExtDtPth, CfgFileName)
CfgFile
}
