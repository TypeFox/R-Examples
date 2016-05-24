
# CONTROL COMPONENT =================================================


# RUNCLASS =====================

#' RunClass 
#' 
#' RunClass is an class contain ReportClass and AnalysisClass objects as separate lists.
#' A RunClass object is the output of running the main function preproviz() and can be
#' studied with either get methods (AnalysisClass objects) or plot methods (ReportClass
#' objects).
#'  
#' @slot reports A list of ReportClass objects
#' @slot analysis A list of AnalysisClass object
#' @return A RunClass object
#' @export

setClass("RunClass", representation(reports="list", analysis="list")) #

#' getconstructeddata RunClass
#' @describeIn getconstructeddata

setMethod("getconstructeddata", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "constructeddata"))
}
)

#' getbasedata RunClass
#' @describeIn getbasedata

setMethod("getbasedata", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "basedata"))
}
)

#' getminmaxconstructeddata RunClass
#' @describeIn getminmaxconstructeddata

setMethod("getminmaxconstructeddata", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "minmaxconstructeddata"))
}
)

#' getcombineddata RunClass
#' @describeIn getcombineddata

setMethod("getcombineddata", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "combineddata"))
}
)

#' getnumericombineddata RunClass
#' @describeIn getnumericombineddata

setMethod("getnumericombineddata", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "numericombineddata"))
}
)

#' getclasslabels RunClass
#' @describeIn getclasslabels

setMethod("getclasslabels", signature(object = "RunClass"), function(object) {
  listcmds <- lapply(object@analysis, function(x) slot(x, "classlabel"))
}
)

# CONTROLCLASS =====================

#' A ControlClass is an class containing SetUpClass objects.  
#'
#' @slot setups A list of SetUpClass objects 

setClass("ControlClass", representation(setups="list"))

#' initializecontrolclassobject
#'
#' initializecontrolclassobject is a constructor function for intializing ControlClass objects. 
#' @param setups (list) Name of SetUpClass objects 
#' @return (ControlClass) object 
#' @export

initializecontrolclassobject <- function(setups){
  if (class(setups)!="list") {stop("Argument to initializecontrolclassobject must be a list.")}
  if (length(setups)==0) {stop("Argument to initializecontrolclassobject must have one or more list elements.")}
  #if (!all(unlist(lapply(setups, function(x) class(eval(as.name(x)))=="SetUpClass")))) {stop("Argument to initializecontrolclassobject must point to SetUpClass object(s).")}
  
  controlclassobject <- new("ControlClass")
  controlclassobject@setups <- setups
  return(controlclassobject)
}

# PARAMETERCLASS =====================

#' ParameterClass 
#' 
#' ParameterClass is a class containing a list of sub class objects (i.e. constructed features, inhereted
#' from BaseClass). A ParameterClass object in a SetUpClass object defines which constructed features are computed 
#' from a DataClass object  
#'
#' @slot parameters A list of sub class objects
#' @export

setClass("ParameterClass", representation(parameters="list"))

#' initializeparameterclassobject
#'
#' initializeparameterclassobject is a constructor function for intializing ParameterClass objects. 
#' @param parameters (list) Name of sub classes 
#' @return (ParameterClass) object 
#' @export

initializeparameterclassobject <- function(parameters){
  if (class(parameters)!="list") {stop("Argument to initializeparameterclassobject must be a list.")}
  if (length(parameters)<3) {stop("Argument to initializeparameterclassobject must have three or more list elements.")}
  if (!all(unlist(lapply(parameters, function(x) extends(x, "BaseClass"))))) {stop("Argument to initializecontrolclassobject must point to the names of defined sub classes inhereted from BaseClass.")}
  
  parameterclassobject <- new("ParameterClass")
  parameterclassobject@parameters <- parameters
  return(parameterclassobject)
}


#' getparameters
#' 
#' getparameters is a method for getting the list of parameters from a ParameterClass object or a Parameter object from SetUpClass object
#' @param object (ParameterClass or SetUpClass)
#' @rdname getparameters
#' @export

setGeneric("getparameters",function(object){standardGeneric("getparameters")}) ## VALUE

#' getparameters ParameterClass
#' @describeIn getparameters

setMethod("getparameters", signature(object = "ParameterClass"), function(object) {
  return(object@parameters)}
)

# DATACLASS =====================

#' DataClass 
#' 
#' DataClass is a class containing data.  
#' DataClass object can be initialized only for a data frame that has a) one class label columns of class 'factor' and
#' b) other columns are of type 'numeric'
#'
#' @slot name (character) name of the setup object
#' @slot basedata (data frame) original data to be visualized 
#' @slot imputedbase (data frame) missing value in original data imputed with Knn imputation 
#' @slot numericdata (data frame) numeric columns of original data
#' @slot imputednumeric (data frame) imputed numeric columns
#' @export

setClass("DataClass", representation(name="character", basedata="data.frame", imputedbase="data.frame", numericdata="data.frame", imputednumeric="data.frame"))

#' initializedataObject 
#' 
#' initializedataobject is a constructor function for initializing aDataClass object from original data frame
#' @param data (data frame)
#' @return (DataClass)
#' @export

initializedataobject <- function(data){

  if(class(data)!="data.frame"){stop("Argument 'data' is not of class data frame.")}
  if(sum(sapply(data, is.factor)==TRUE)!=1) {stop("Argument to initializedataobject must have one and only one factor column.")}
  if(sum(sapply(data, is.numeric)==TRUE)!=ncol(data)-1) {stop("Argument initializedataobject must have only numeric columns and one factor column.")}
  
  classlabels <- data[sapply(data, is.factor)][,1]
  if (any(is.na(classlabels))==TRUE) {stop("Factor column in data can not have missing values.")}
  
  dataobject <- new("DataClass", basedata=data)
  dataobject@name <- as.character(substitute(data))
  suppressWarnings(dataobject@imputedbase <- DMwR::knnImputation(data, k=5))
  dataobject@imputednumeric <- dataobject@imputedbase[sapply(dataobject@imputedbase,is.numeric)]
  dataobject@numericdata <- data[sapply(data,is.numeric)]

  return(dataobject)

}

# SETUPCLASS =====================

#' SetUpClass 
#' 
#' SetUpClass is an class containing a DataClass object, a ParameterClass object and the name of the SetUpClass object  
#'
#' @slot data (DataClass) 
#' @slot parameters (ParameterClass) 
#' @slot objectname (character)
#' @export

setClass("SetUpClass", representation(objectname="character", data="DataClass", parameters="ParameterClass"))

#' initializesetupclassobject 
#' 
#' initializesetupclassobject is a constructor function for initializing a SetUpClass object
#' @param objectname (character) Name of the setup
#' @param parameterobject (ParameterClass)
#' @param dataobject (DataClass)
#' @export

initializesetupclassobject <- function(objectname, parameterobject, dataobject){
  if (class(objectname)!="character") {stop("Argument objectname to initializesetupclassobject must be of class character.")}
  if (class(parameterobject)!="ParameterClass") {stop("Argument parameterobject to initializesetupclassobject must be of class ParameterClass.")}
  if (class(dataobject)!="DataClass") {stop("Argument dataobject to initializesetupclassobject must be of class DataClass.")}
  
  setupclassobject <- new("SetUpClass")
  setupclassobject@objectname <- objectname
  setupclassobject@parameters <- parameterobject
  setupclassobject@data <- dataobject
  return (setupclassobject)
}


#' getparameters SetUpClass
#' @describeIn getparameters

setMethod("getparameters", signature(object = "SetUpClass"), function(object) {
  return(object@parameters)}
)

#setGeneric("saveSetUp",function(object){standardGeneric("saveSetUp")})
#
#setMethod("saveSetUp", signature(object = "SetUpClass"), function(object) {
#  fn <- tempfile()
#  save(object, ascii=FALSE, file=fn)
#}
#  )




###





