##' @import utils
##' @export
`CallMultiSV` <- function(CfgFile, LgLim, SVSize, MtlSVOut) UseMethod("CallMultiSV")

.onAttach <- 
		function(lib, pkg, ...)
{
	MltSvDiscp <- packageDescription(pkg)
	MltSvVer <- MltSvDiscp$Version
	MltSvDate <- MltSvDiscp$Date
	MltSvNm <- MltSvDiscp$Package
	MltSvTl <- MltSvDiscp$Title
	MltSvAthr <- MltSvDiscp$Author
	MltSvMntr <- MltSvDiscp$Maintainer
	packageStartupMessage(paste("\n", MltSvNm, ": ", MltSvTl, sep = ""))
	packageStartupMessage(paste("Version ", MltSvVer, " (", MltSvDate, ") ", sep = ""))
	packageStartupMessage(paste("Authors: ", MltSvAthr, sep = ""))
	packageStartupMessage(paste("Maintainer: ", MltSvMntr, "\n", sep = ""))
	packageStartupMessage('Use citation("MultiSV") to know how to cite our work.\n')
}
##' @export
`ComputeBinCounts` <- function(RDBinSampleFile, RDBinChrSizeFile, RDBinSampleName , 
					 RDBinWindowSize
						, OutFolder) UseMethod("ComputeBinCounts") 

##' @export
`CmptMltPvl` <- function (MultiData,CfgFile) UseMethod("CmptMltPvl") 

#' @import reshape
##' @export

`IdfMltSV` <- function(MultiData,LgLim,SVSize) UseMethod("IdfMltSV") 

#' @import reshape
##' @export

`WriteMultiSV` <- function(MultiData,File) UseMethod("WriteMultiSV") 

##' @import nlme
##' @export

`PrbMlt` <- function(MultiData) UseMethod("PrbMlt")

##' @export

###############################################################################
`gtExEs` <- function(proc) UseMethod("gtExEs")

##' @export
`Bin2MultiSV` <- function(CfgFile) UseMethod("Bin2MultiSV")

##' @export
`PrepMultiDtLgMn` <- function (MultiData) UseMethod("PrepMultiDtLgMn")

##' @export
`PrepMultiDt` <- function (MultiData) UseMethod("PrepMultiDt")


##' @export
`ProcMutilDt` <- function (MultiData) UseMethod("ProcMutilDt")


##' @export
`MultiSVExample` <- function (MultiSVData) UseMethod("MultiSVExample")




