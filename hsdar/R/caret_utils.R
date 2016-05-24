  
.getAttrParameters <- function()  
  c("response", "predictor")
  
.getCaretCompatibleClasses <- function()
  c("Speclib", "Nri")
 
setClassUnion(".CaretHyperspectral", .getCaretCompatibleClasses())

setOldClass("preProcess")

setClass(".preProcessHyperspectral",
         representation(
           preProcess = "preProcess"
         ),
         prototype(
           preProcess = list()
         )
)

setMethod("show", ".preProcessHyperspectral", function(object) show(object@preProcess))