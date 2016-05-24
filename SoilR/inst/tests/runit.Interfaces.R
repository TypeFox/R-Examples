#!/usr/bin/Rscript
# vim:set ff=unix expandtab ts=2 sw=2:
test.TimeMapInterface <- function(){
  # 
  classNames=c(
    "TimeMap",
    "TransportDecompositionOperator",
    "BoundLinDecompOp",
    "ConstLinDecompOp",
    "BoundInFlux",
    "BoundFc")
  methodNames=c(
    "getTimeRange",
    "getFunctionDefinition"
  )
  m=matrix(nrow=length(classNames),ncol=length(methodNames))
  colnames(m)<-methodNames
  rownames(m)<-classNames
  for (f in  methodNames){
    for (c in classNames) {
      m[c,f]=hasMethod(f=f,signature=c)
    }
  }
  val <- all(m)
  if (!val){print(m)}
  checkTrue(val)

}
