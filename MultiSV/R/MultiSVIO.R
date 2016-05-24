##' @export
`CallMultiSV.default` <- function(CfgFile,LgLim,SVSize,MtlSVOut){
MultiData <- Bin2MultiSV(CfgFile)
MultiData <- PrepMultiDt(MultiData)
MultiData <- PrepMultiDtLgMn(MultiData)
MultiData <- ProcMutilDt(MultiData)
WriteMultiSV(IdfMltSV(MultiData,LgLim,SVSize),MtlSVOut)
}
