setClassUnion('dataframeOrNULL', c('data.frame', 'NULL'))
setClassUnion('numericOrNULL', c('numeric', 'NULL'))
setClassUnion('characterOrNULL', c('character', 'NULL'))
setClassUnion('logicalOrNULL', c('logical', 'NULL'))
setClassUnion('matrixOrNULL', c('matrix', 'NULL'))
setClassUnion('listOrNULL', c('list', 'NULL'))

#' @useDynLib simPop
#' @import doParallel
#' @import foreach
#' @import nnet
#' @import data.table
#' @import e1071
#' @import parallel
#' @import vcd
#' @import methods
#' @import Rcpp
#' @import party
#' @importFrom lattice bwplot panel.bwplot packet.number panel.points panel.refline panel.xyplot
#' @importFrom laeken calibVars
#' @importFrom MASS ginv
#' @importFrom colorspace heat_hcl
#' @importFrom VIM hotdeck
#' @importFrom graphics par
#' @importFrom stats as.formula chisq.test coef cor cov dlogis formula lm mad quantile rexp
#' @importFrom utils read.csv2 tail head
#' @importFrom stats median model.matrix optim plogis ppoints predict rnorm runif uniroot var weighted.mean glm poisson
# removed  stats quantilerexp
#' @importFrom plyr revalue
NULL



#' Class \code{"dataObj"}
#'
#' Objects of this class contain information on a population or survey.
#'
#'
#' @name dataObj-class
#' @aliases dataObj-class show,dataObj-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("dataObj", ...)} but are usually automatically created when using
#' \code{\link{simStructure}}.
#' @author Bernhard Meindl and Matthias Templ
#' @seealso \code{\linkS4class{simPopObj}}
#' @keywords classes
#' @examples
#' showClass("dataObj")
#'
#' ## show method, generate an object of class dataObj first
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", weight="rb050", strata="db040")
#' ## shows some basic information:
#' inp
#'
NULL

#' Extract and modify variables from population or sample data stored in an
#' object of class \code{\link{simPopObj-class}}.
#'
#' Using \code{\link{samp}} \code{\link{samp<-}} it is possible to extract or
#' rather modify variables of the sample data within slot \code{data} in slot
#' \code{sample} of the \code{\link{simPopObj-class}}-object. Using
#' \code{\link{pop}} \code{\link{pop<-}} it is possible to extract or rather
#' modify variables of the synthetic population within in slot \code{data} in
#' slot \code{sample} of the \code{\link{simPopObj-class}}-object.
#'
#'
#' @name get_set-methods
#' @aliases pop pop<- pop,simPopObj-method pop<-,simPopObj-method
#' samp samp<- samp,simPopObj-method samp<-,simPopObj-method
#' popData popData,simPopObj-method
#' sampleData,simPopObj-method
#' popObj popObj<- popObj<-,simPopObj,dataObj-method popObj,simPopObj-method
#' sampleData
#' sampleObj sampleObj<- sampleObj,simPopObj-method sampleObj<-,simPopObj,dataObj-method
#' tableObj tableObj,simPopObj-method
#' @docType methods
#' @param obj An object of class \code{\link{simPopObj-class}}
#' @param var variable name or index for the variable in slot 'samp' of object
#' with the slot name to be accessed. If \code{NULL}, the entire dataset
#' (sample or population) is returned.
#' @param value Content replacing whatever the variable in slot \code{var} in
#' \code{obj} currently holds.
#' @return Returns an object of class \code{\link{simPopObj-class}} with the
#' appropriate replacement.
#' @author Bernhard Meindl
#' @seealso \code{\link{simPopObj-class}},\code{\link{pop}},
#' \code{\link{pop<-}}, \code{\link{samp<-}}, \code{\link{manageSimPopObj}}
#' @keywords manip methods
#' @examples
#'
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040",
#' weight="db090")
#' simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#'
#' ## get/set variables in sample-object of simPopObj
#' head(samp(simPopObj, var="age"))
#' samp(simPopObj, var="newVar") <- 1
#' head(samp(simPopObj, var="newVar"))
#' ## deleting is also possible
#' samp(simPopObj, var="newvar") <- NULL
#' head(samp(simPopObj, var="newvar"))
#' ## extract multiple variables
#' head(samp(simPopObj, var=c("db030","db040")))
#'
#' ## get/set variables in pop-object of simPopObj
#' head(pop(simPopObj, var="age"))
#' pop(simPopObj, var="newVar") <- 1
#' head(pop(simPopObj, var="newVar"))
#' ## deleting is also possible
#' pop(simPopObj, var="newvar") <- NULL
#' head(pop(simPopObj, var="newvar"))
#' ## extract multiple variables
#' head(pop(simPopObj, var=c("db030","db040")))
#'
NULL



#' Class \code{"simPopObj"}
#'
#' An object that is used throughout the package containing information on the
#' sample (in slot \code{sample}), the population (slot \code{pop}) and
#' optionally some margins in form of a table (slot \code{table}).
#'
#'
#' @name simPopObj-class
#' @aliases simPopObj-class show,simPopObj-method
#' @docType class
#' @section Objects from the Class: Objects are automatically created in
#' function \code{\link{simStructure}}.
#' @author Bernhard Meindl and Matthias Templ
#' @seealso \code{\linkS4class{dataObj}}
#' @keywords classes
#' @examples
#'
#' showClass("simPopObj")
#'
#' ## show method: generate an object of class simPop first
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' class(eusilcP)
#' ## shows some basic information:
#' eusilcP
#'
NULL


setClass(
  Class='dataObj',
  representation=representation(
    data='dataframeOrNULL',
    hhid='characterOrNULL',
    pid='characterOrNULL',
    hhsize='characterOrNULL',
    weight='characterOrNULL',
    strata='characterOrNULL',
    ispopulation='logicalOrNULL'
  ),
  prototype=prototype(
    data=NULL,
    hhid=NULL,
    pid=NULL,
    hhsize=NULL,
    weight=NULL,
    strata=NULL,
    ispopulation=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)

setClassUnion('dataObjOrNULL', c('dataObj', 'NULL'))

setClass(
  Class='simPopObj',
  representation=representation(
    sample='dataObjOrNULL',
    table='dataframeOrNULL',
    pop='dataObjOrNULL',
    basicHHvars='characterOrNULL'
  ),
  prototype=prototype(
    sample=NULL,
    table=NULL,
    pop=NULL,
    basicHHvars=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)

#' get and set variables from population or sample data stored in an object of
#' class \code{\linkS4class{simPopObj}}.
#'
#' This functions allows to get or set variables in slots \code{pop} and
#' \code{sample} of \code{\linkS4class{simPopObj}}-objects. This is a utility
#' function that is useful for writing custom wrapper functions.
#'
#' @name manageSimPopObj
#' @param x an object of class \code{\linkS4class{simPopObj}}.
#' @param var character vector of length 1; variable name that should be set or
#' extracted.
#' @param sample a logical indicating whether \code{var} should be
#' extracted/set from slot 'sample' (TRUE) or slot 'pop' (FALSE).
#' @param set logical; if TRUE, argument 'values' is set to either the sample
#' or population data stored in 'x', depending on argument 'sample'. If FALSE,
#' the desired variable given by 'var' is returned from either the sample or
#' the pop slot of 'x'.
#' @param values vector; if 'set' is TRUE, then this vector is used to update
#' the variable of sample or population data depending of choice of argument
#' 'sample'.
#' @return An object of class \code{\linkS4class{simPopObj}} (if 'set' is TRUE)
#' or a vector (if 'set' is FALSE).
#' @export
#' @author Bernhard Meindl and Matthias Templ
#' @keywords manip
#' @examples
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040",
#'   weight="db090")
#' simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#'
#' (manageSimPopObj(simPopObj, var="age", sample=FALSE, set=FALSE))
#' (manageSimPopObj(simPopObj, var="age", sample=TRUE, set=FALSE))
manageSimPopObj <- function(x, var, sample=FALSE, set=FALSE, values=NULL) {
  if ( class(x) != "simPopObj" ) {
    stop("wrong input of argument 'x' (needs to be of class 'simPopObj')!\n")
  }
  if ( length(var) != 1 ) {
    stop("only one variable can be used at a time!\n")
  }
  if ( set==FALSE ) {
    if ( sample ) {
      return(invisible(samp(x, var=var)))
    } else {
      return(invisible(pop(x, var=var)))
    }
  }
  if ( set == TRUE ) {
    if ( is.null(values) ) {
      stop("you need to provide values!\n")
    }
    if ( sample ) {
      samp(x, var=var) <- values
    } else {
      pop(x, var=var) <- values
    }
    return(invisible(x))
  }
}

#' @export
setGeneric("samp", function(obj, var=NULL) {
  standardGeneric("samp")
})
setMethod("samp", "simPopObj", function(obj, var=NULL) {
  sampData <- as.data.table(obj@sample@data)
  cn <- names(sampData)
  if ( is.null(var) ) {
    return(sampData)
  }
  if ( is.numeric(var) ) {
    varI <- match(var, 1:ncol(sampData))
    vars <- cn[which(!is.na(varI))]
  }
  if ( is.character(var) ) {
    varI <- match(var, names(sampData))
    var <- var[which(!is.na(varI))]
  }
  if ( all(is.na(varI)) ) {
    warning("No specified variable was found in the sample! Check input 'var'!\n")
    return(NULL)
  }
  if ( any(is.na(varI)) ) {
    ww <- "The following variables/indices were not found in the sample:\n"
    ww <- paste0(ww, paste(var[is.na(varI)], collapse=" | "))
    warning(ww)
  }
  
  if ( length(var)==1 ) {
    return(sampData[[var]])
  } else {
    return(sampData[,var,with=F])
  }
})

#' @export
setGeneric("samp<-", function(obj, var, value) {
  standardGeneric("samp<-")
})
#' @export
setReplaceMethod("samp", "simPopObj", function(obj, var, value) {
  if ( length(var) != 1) {
    stop("we can only set one variable!\n")
  }
  obj@sample@data[[var]] <- value
  validObject(obj)
  obj
})

#' @export
setGeneric("pop", function(obj, var=NULL) {
  standardGeneric("pop")
})

#' @export
setMethod("pop", "simPopObj", function(obj, var=NULL) {
  popData <- as.data.table(obj@pop@data)
  cn <- names(popData)
  if ( is.null(var) ) {
    return(popData)
  }
  if ( is.numeric(var) ) {
    varI <- match(var, 1:ncol(popData))
    vars <- cn[which(!is.na(varI))]
  }
  if ( is.character(var) ) {
    varI <- match(var, names(popData))
    var <- var[which(!is.na(varI))]
  }
  if ( all(is.na(varI)) ) {
    warning("No specified variable was found in the synthetic population! Check input 'var'!\n")
    return(NULL)
  }
  if ( any(is.na(varI)) ) {
    ww <- "The following variables/indices were not found in the synthetic population:\n"
    ww <- paste0(ww, paste(var[is.na(varI)], collapse=" | "))
    warning(ww)
  }
  if ( length(var)==1 ) {
    return(popData[[var]])
  } else {
    return(popData[,var,with=F])
  }
})

#' @export
setGeneric("pop<-", function(obj, var, value) {
  standardGeneric("pop<-")
})
#' @export
setReplaceMethod("pop", "simPopObj", function(obj, var, value) {
  if ( length(var) != 1) {
    stop("we can only set one variable!\n")
  }
  obj@pop@data[[var]] <- value
  validObject(obj)
  obj
})

#' @export
setGeneric("sampleObj", function(object) standardGeneric("sampleObj"))
#' @export
setGeneric("sampleObj<-", function(object, value) standardGeneric("sampleObj<-"))
#' @export
setMethod("sampleObj", signature="simPopObj", definition=function(object) {
  invisible(object@sample)
})
#' @export
setReplaceMethod("sampleObj", signature = c("simPopObj", "dataObj"), definition = function(object, value) {
  object@sample <- value
  validObject(object)
  invisible(object)
})

#' @export
setGeneric("sampleData", function(object) standardGeneric("sampleData"))
#' @export
setMethod("sampleData", signature="simPopObj", definition=function(object) {
  invisible(sampleObj(object)@data)
})

#' @export
setGeneric("popObj", function(object) standardGeneric("popObj"))
#' @export
setGeneric("popObj<-", function(object, value) standardGeneric("popObj<-"))
#' @export
setMethod("popObj", signature="simPopObj", definition=function(object) {
  invisible(object@pop)
})
#' @export
setReplaceMethod("popObj", signature = c("simPopObj", "dataObj"), definition = function(object, value) {
  object@pop <- value
  validObject(object)
  invisible(object)
})

#' @export
setGeneric("popData", function(object) standardGeneric("popData"))
#' @export
setMethod("popData", signature="simPopObj", definition=function(object) {
  invisible(popObj(object)@data)
})

#' @export
setGeneric("tableObj", function(object) standardGeneric("tableObj"))
#' @export
setMethod("tableObj", signature="simPopObj", definition=function(object) {
  invisible(object@table)
})
