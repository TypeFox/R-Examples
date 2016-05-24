#' query \code{safeObj}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{safeObj}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item dimInfo: get infos on dimensional variables that formed the base of the protected data
#' \item elapsedTime: get elapsed time of the protection procedure
#' \item finalData: return final data object
#' \item nrNonDuplicatedCells: total number of cells that are duplicates
#' \item nrPrimSupps: total number of primary suppressed cells
#' \item nrSecondSupps: total number of secondary cell suppressions
#' \item nrPublishableCells: total number of cells that can be published
#' \item suppMethod: suppression method that has been used
#' \item cellInfo: extract information about a specific cell
#' \item cellID: calculate ID of a specific cell defined by level-codes and variable names
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type matches 'dimInfo', 'elapsedTime', 'finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells' or 'suppMethod': input is not used (empty list)
#' \item type matches 'cellInfo' or 'cellID': input is a list of length 3
#' \itemize{
#' \item first element: character vector specifying variable names that need to exist in slot 'dimInfo' of \code{object}
#' \item second element: character vector specifying codes for each variable that define a specific table cell
#' \item third element: logical vector of length 1 with TRUE setting verbosity and FALSE to turn verbose output off}
#'
#' @return information from \code{object} depending on \code{type}
#' \itemize{
#' \item an object of class \code{dimInfo} (or NULL) if type matches 'dimInfo'
#' \item a numeric vector if type matches 'elapsedTime', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells' or 'cellID'
#' \item a character vector if type matches 'suppMethod'
#' \item a data.frame if type matches 'finalData'
#' \item a list if type matches 'cellInfo' containing the following elements:
#' \itemize{
#' \item element 'cellID': numeric vector of length 1 specifying the index of the cell of interest
#' \item element 'data': row of slot 'finalData' with the row being defined by the calculated \code{cellID}
#' \item element 'primSupp': logical vector of length 1 being TRUE if cell is a primary suppressed cell
#' \item element 'secondSupps': logical vector of length 1 being TRUE if cell is a secondary suppressed cell }
#' }
#'
#' @export
#' @docType methods
#' @rdname get.safeObj-method
#' @aliases get.safeObj,safeObj,character,list-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("get.safeObj", function(object, type, input) { standardGeneric("get.safeObj")})

setGeneric("g_dimInfo", function(object) { standardGeneric("g_dimInfo") })
setGeneric("g_elapsedTime", function(object) { standardGeneric("g_elapsedTime") })
setGeneric("g_finalData", function(object) { standardGeneric("g_finalData") })
setGeneric("g_nrNonDuplicatedCells", function(object) { standardGeneric("g_nrNonDuplicatedCells") })
setGeneric("g_nrPrimSupps", function(object) { standardGeneric("g_nrPrimSupps") })
setGeneric("g_nrSecondSupps", function(object) { standardGeneric("g_nrSecondSupps") })
setGeneric("g_nrPublishableCells", function(object) { standardGeneric("g_nrPublishableCells") })
setGeneric("g_suppMethod", function(object) { standardGeneric("g_suppMethod") })
setGeneric("g_getCellInfo", function(object, input) { standardGeneric("g_getCellInfo") })
setGeneric("g_getCellID", function(object, input) { standardGeneric("g_getCellID") })
