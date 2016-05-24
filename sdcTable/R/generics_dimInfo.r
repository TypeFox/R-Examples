#' query \code{dimInfo}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dataObj}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item strInfo: info on how many digits in the default codes ach dimensional variable allocates
#' \item dimInfo: a list object with each slot containing an object of class \code{dimVar}
#' \item varName: variable names
#' \item strID: character vector of ID's defining table cells
#' \item posIndex vector showing the index of the elements of \code{dimInfo} in the underlying data
#'
#' @return information from objects of class \code{dimInfo} depending on argument \code{type}
#' \itemize{
#' \item a list (or NULL) if argument \code{type} matches 'strInfo', 'dimInfo'
#' \item numeric vector (or NULL) if argument \code{type} matches 'posIndex'
#' \item character vector (or NULL) if argument \code{type} matches 'varName' or 'strID'}
#'
#' @export
#' @docType methods
#' @rdname get.dimInfo-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("get.dimInfo", function(object, type) { standardGeneric("get.dimInfo") })

#' modify \code{dimInfo}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dimInfo}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item strID: set slot 'strID' of argument \code{object}
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type==strID: a character vector containing ID's
#'
#' @return an object of class \code{dimInfo}
#'
#' @export
#' @docType methods
#' @rdname set.dimInfo-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("set.dimInfo", function(object, type, input) { standardGeneric("set.dimInfo") })

# get-methods
setGeneric("g_str_info", function(object) { standardGeneric("g_str_info") })
setGeneric("g_dim_info", function(object) { standardGeneric("g_dim_info") })
setGeneric("g_varname", function(object) { standardGeneric("g_varname") })
setGeneric("g_str_id", function(object) { standardGeneric("g_str_id") })
setGeneric("g_pos_index", function(object) { standardGeneric("g_pos_index") })

# set-methods
setGeneric("s_str_id<-", function(object, value) { standardGeneric("s_str_id<-") })

