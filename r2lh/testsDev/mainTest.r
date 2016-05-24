
if(getwd()=="C:/Documents and Settings/Christophe/Mes documents"){
    setwd("C:/Documents and Settings/Christophe/Mes documents/Recherche/R et SAS/r2lBiv/trunk/Xolh/testsDev")
}else{}
if(getwd()=="C:/Documents and Settings/Administrator/My Documents"){
    setwd("C:/Documents and Settings/Administrator/My Documents/Recherche/R et SAS/yeap/trunk/Xolh/testsDev")
}else{}

library(codetools)

cleanProg <- function(realResult,theoResult="",result=TRUE,tolerance=0){
  functionNames <- strsplit(deparse(substitute(realResult)),"\\(")[[1]][1]
  if(identical(theoResult,"")==FALSE){
    if( isTRUE(all.equal( realResult , theoResult ))!=result ){
      cat("WARNING(PreTest2) in    ",functionNames,":",deparse(substitute(realResult)), " == ",theoResult," is not ",result,"\a\n\a")
    }
  }else{}
  if(length(findGlobals(get(functionNames),FALSE)$variables)  > tolerance){
    cat("WARNIGS(detectGlobal) in ",functionNames,": These are the globals:",findGlobals(get(functionNames),FALSE)$variables,"\a\n")
    stop()
  }else{}
}


source("testDefineVariable.r")
source("testFunctions.r")

source("testRtlu.r")

source("testLogical.r")
source("testFactor.r")
source("testOrdered.r")
source("testDiscrete.r")
source("testContinuous.r")

#source("testDisplay.r")
source("testRtlb.r")
