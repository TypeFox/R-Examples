#' query \code{dataObj}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dataObj}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item rawData: raw input data
#' \item dimVarInd: indices of dimensional variables
#' \item freqVarInd: index of frequency variable
#' \item numVarInd: indices of numerical variables
#' \item weightVarInd index of weight variable
#' \item sampWeightInd index of variable holding sampling weights
#' \item isMicroData does \code{object} consist of microdata?
#' \item numVarNames variable names of numerical variables
#' \item freqVarName variable name of frequency variable
#' \item varName variable names of dimensional variables
#'
#' @return information from objects of class \code{dataObj} depending on argument \code{type}
#' \itemize{
#' \item a list if argument \code{type} matches 'rawData'
#' \item numeric vector if argument \code{type} matches 'dimVarInd', 'freqVarInd', 'numVarInd', 'weightVarInd' or 'sampWeightInd'
#' \item character vector if argument \code{type} matches 'numVarNames', 'freqVarName' or 'varName'
#' \item logical vector of length 1 if argument \code{type} matches 'isMicroData'}
#'
#' @export
#' @docType methods
#' @rdname get.dataObj-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("get.dataObj", function(object, type) { standardGeneric("get.dataObj")})

#' modify \code{dataObj}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dataObj}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item rawData: set slot 'rawData' of argument \code{object}
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type==rawData: a list containing raw data
#'
#' @return an object of class \code{dataObj}
#'
#' @export
#' @docType methods
#' @rdname set.dataObj-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("set.dataObj", function(object, type, input) { standardGeneric("set.dataObj") })

#' initialize \code{dataObj}-objects
#'
#' @param input a list with element described below:}
#'
#' \itemize{
#' \item element 'inputData': a list object holding data
#' \item element 'dimVarInd': index (within \code{inputData}) of variables that define the table to protect
#' \item element 'freqVarInd': index (within \code{inputData}) of variable holding frequencies
#' \item element 'numVarInd' index (within \code{inputData}) of numerical variables (or NULL)
#' \item element 'weightInd': index (within \code{inputData}) of variable holding weights (or NULL)
#' \item element 'sampWeightInd': index (within \code{inputData}) of variable holding sampling weights (or NULL)
#' \item element 'isMicroData': logical vector of length 1
#'
#' @return an object of class \code{dataObj}
#'
#' @export
#' @docType methods
#' @rdname init.dataObj-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("init.dataObj", function(input) { standardGeneric("init.dataObj") })

# get-methods
setGeneric("g_raw_data", function(object) { standardGeneric("g_raw_data") })
setGeneric("g_dimvar_ind", function(object) { standardGeneric("g_dimvar_ind") })
setGeneric("g_freqvar_ind", function(object) { standardGeneric("g_freqvar_ind") })
setGeneric("g_numvar_ind", function(object) { standardGeneric("g_numvar_ind") })
setGeneric("g_weightvar_ind", function(object) { standardGeneric("g_weightvar_ind") })
setGeneric("g_sampweight_ind", function(object) { standardGeneric("g_sampweight_ind") })
setGeneric("g_is_microdata", function(object) { standardGeneric("g_is_microdata") })
setGeneric("g_numvar_names", function(object) { standardGeneric("g_numvar_names") })
setGeneric("g_freqvar_name", function(object) { standardGeneric("g_freqvar_name") })
setGeneric("g_var_name", function(object) { standardGeneric("g_var_name") })

# set-methods
setGeneric("s_raw_data<-", function(object, value) { standardGeneric("s_raw_data<-") })

