###
### AllClasses.R
###

#require(methods)

##=============================================================================
setClass("msPriorSpec",
         representation(priorType="character",
                        priorDistr="character",
                        priorPars="vector"),
         prototype(priorPars=NA))


##=============================================================================
setClass("msfit", representation("list"), prototype = prototype(elementType = "list"), contains="list")
