##' @export
`MultiSVExample.default` <- function(MultiSVData){
  LgLim = 0.6
  SVSize = 5000
  MultiData <- PrepMultiDt(MultiSVData)
  MultiData <- PrepMultiDtLgMn(MultiData)
  MultiData <- ProcMutilDt(MultiData)
  IdfMltSV(MultiData,LgLim,SVSize)
}
