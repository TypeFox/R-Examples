
# BASECLASS DEFINITION AND GENERIC FUNCTION COMPUTEVALUE ======================

#' @include 00Utils.R
NULL

#' BaseClass
#' 
#' A BaseClass is an abstract class. Sub classes are inherited from BaseClass. 
#'
#' @slot objectname (character) name of the object
#' @slot valuevector (numeric) constructed feature vector
#' @slot isvalid (logical) result of object validation: checksubclassobjectvalidity
#' @slot preimpute (logical) whether valuevector should be compute before missing value imputation
#' @export

setClass("BaseClass", representation(objectname="character", valuevector="numeric", isvalid="logical", preimpute="logical"),
         prototype(isvalid=FALSE))

#' computeValue
#' 
#' computevalue is a BaseClass generic function for computing constructed feature vectors.
#' @param object A sub class object inheriting BaseClass 
#' @param dataobject A DataClass object
#' @return (numeric) feature vector
#' @export

setGeneric("computeValue", function(object, dataobject) {
  standardGeneric("computeValue")
})

# SUB CLASS DEFINITION AND METHOD DEFINITION ====================

#' constructfeature
#'
#' constructfeature is a convenience function for inheriting a user-defined sub class from BaseClass and defining computevalue methods for them.
#' @param classname (character) Name of the inherited class 
#' @param operation (expression) Feature construction operation . The expression is evaluated by computevalue method. 
#' @param mode (character) Mode of data to be used in construction . Defaults to "all", option "numeric" for numeric data without class labels. 
#' @param impute (logical) Impute whether construction is done before missing value imputation . Defaults to "FALSE"
#' @return (BaseClass) A sub class inhereted from BaseClass
#' @return A computevalue method for the sub class
#' @export

constructfeature <- function(classname, operation, mode="all", impute=FALSE){

  setClass(classname, contains="BaseClass", where=topenv(parent.frame()), prototype=prototype(preimpute=impute))

  if (mode=="all" & impute=="TRUE") {functionexpression <- gsub("data", "dataobject@basedata", operation)}
  if (mode=="numeric" & impute=="TRUE") {functionexpression <- gsub("data", "dataobject@numericdata", operation)}
  if (mode=="all" & impute=="FALSE") {functionexpression <- gsub("data", "dataobject@imputedbase", operation)}
  if (mode=="numeric" & impute=="FALSE") {functionexpression <- gsub("data", "dataobject@imputednumeric", operation)}

  setMethod("computeValue", where=topenv(parent.frame()), signature(object = classname), function(object, dataobject) {
    #print(paste("Computing feature vector:", object@objectname, " from data: ", dataobject@name, sep=""))
    temp <- object ## BUG: This is completely redundant
    eval(parse(text=functionexpression))
  })

}

# VALIDITY ===============================

#' checksubclassobjectvalidity
#' 
#' checksubclassobjectvalidity is a validation method for BaseClass objects.
#' Valid constructed feature vectors must be finite, have variance and be of same length as rows in data.
#' 
#' @param object A user-defined or default sub class object inhereted from BaseClass
#' @param dataobject An object of DataClass 
#' @return (numeric) 1 for valid objects
#' @return (numeric) 0 for not valid object

checksubclassobjectvalidity <- function(object, dataobject) {
  
  errors <- character()

  isfinite_valuevector <- all(is.finite(object@valuevector))
  if (isfinite_valuevector == FALSE) {
    msg <- paste("Computation has elements that are not finite: Inf, NaN, NA, Should be finite.", sep = "")
    errors <- c(errors, msg)
  }
  
  isnearzerovar_valuevector <- length(caret::nearZeroVar(object@valuevector))
  if (isnearzerovar_valuevector != 0) {
    msg <- paste(class(object)[1], ": Constructed feature vector has zero or near zero variance. Object is not valid in this data.", sep = "")
    errors <- c(errors, msg)
  }

  length_valuevector <- length(object@valuevector)==nrow(dataobject@basedata)
  if (length_valuevector != TRUE) {
    msg <- paste("Computation result lenght differs from nrow(basedata)", sep = "")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) {
    # object is valid
    return(1)}
  else {
      return(0)
    }
}

# SUB CLASS OBJECT INITIALIZATION ============================

#' initializesubclassobject
#' 
#' initializesubclassobject is a constructor function for constructing objects for classes inhereted from BaseClass
#' @param classname (character) A name of class of an object is initialized from (e.g. "MissingValueShare").
#' @param dataobject A DataClass object that is used by computeValue method for the sub class object.
#' @export

initializesubclassobject <- function(classname, dataobject){
  y <- new(classname)
  y@objectname <- gsub("[[:punct:]]", "", as.character(deparse(classname)))
  y@valuevector <- computeValue(y, dataobject)
  if (checksubclassobjectvalidity(y, dataobject)==1) {y@isvalid <- TRUE}
  return(y)
}

#' getinitializedsubclassobjects
#' 
#' getinitializedsubclassobjects is a constructor function for constructing sub class objects
#' @param dataobject A DataClass object that is used by computeValue method for the sub class object.
#' @param parameterobject A ParameterClass object that is composed of sub class names

getinitializedsubclassobjects <- function(dataobject, parameterobject)
{
  cls <- parameterobject@parameters
  cls2 <- lapply(cls, function(x) initializesubclassobject(x, dataobject))
}








