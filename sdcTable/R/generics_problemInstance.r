#' query \code{problemInstance}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{problemInstance}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item strID: vector of unique IDs for each table cell
#' \item nrVars: total number of table cells
#' \item freq: vector of frequencies
#' \item w: a vector of weights used in the linear problem (or NULL)
#' \item numVars: a list containing numeric vectors containing values for numerical variables for each table cell (or NULL)
#' \item sdcStatus: a vector containing the suppression state for each cell (possible values are 'u': primary suppression, 'x': secondary suppression, 'z': forced for publication, 's': publishable cell, 'w': dummy cells that are considered only when applying the simple greedy heuristic to protect the table)
#' \item lb: lower bound assumed to be known by attackers for each table cell
#' \item ub: upper bound assumed to be known by attackers for each table cell
#' \item LPL: lower protection level required to protect table cells
#' \item UPL: upper protection level required to protect table cells
#' \item SPL: sliding protection level required to protect table cells
#' \item primSupps: vector of indices of primary sensitive cells
#' \item secondSupps: vector of indices of secondary suppressed cells
#' \item forcedCells: vector of indices of cells that must not be suppressed
#' \item hasPrimSupps: shows if \code{object} has primary suppressions or not
#' \item hasSecondSupps: shows if \code{object} has secondary suppressions or not
#' \item hasForcedCells: shows if \code{object} has cells that must not be suppressed
#' \item weight: gives weight that is used the suppression procedures
#' \item suppPattern: gives the current suppression pattern
#'
#' @return information from objects of class \code{dataObj} depending on argument \code{type}
#' \itemize{
#' \item a list (or NULL) if argument \code{type} matches 'numVars'
#' \item numeric vector if argument \code{type} matches 'freq', 'lb', 'ub', 'LPL', 'UPL', 'SPL', 'weight', 'suppPattern'
#' \item numeric vector (or NULL) if argument \code{type} matches 'w', 'primSupps', 'secondSupps', 'forcedCells'
#' \item character vector if argument \code{type} matches 'strID', 'sdcStatus', ''
#' \item logical vector of length 1 if argument \code{type} matches 'hasPrimSupps', 'hasSecondSupps', 'hasForcedCells'
#' \item numerical vector of length 1 if argument \code{type} matches 'nrVars'
#' }
#'
#' @export
#' @docType methods
#' @rdname get.problemInstance-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('get.problemInstance', function(object, type) {standardGeneric('get.problemInstance')})

#' modify \code{problemInstance}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{problemInstance}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item lb: set assumed to be known lower bounds
#' \item ub: set assumed to be upper lower bounds
#' \item LPL: set lower protection levels
#' \item UPL: set upper protection levels
#' \item SPL: set sliding protection levels
#' \item sdcStatus: change anonymization status
#' @param input a list with elements 'indices' and 'values'.}
#'
#' \itemize{
#' \item element 'indices': numeric vector defining the indices of the cells that should be modified
#' \item element 'values': numeric vector whose values are going to replace current values for cells defined by 'indices' depending on argument \code{type}
#'
#' @return an object of class \code{problemInstance}
#'
#' @export
#' @docType methods
#' @rdname set.problemInstance-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('set.problemInstance', function(object, type, input) {standardGeneric('set.problemInstance')})

#' perform calculations on \code{problemInstance}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{problemInstance}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item makeMasterProblem: create the master problem that is the core of the secondary cell suppression problem
#' \item isProtectedSolution: check if a solution violates any required (upper|lower|sliding) protection levels
#' @param input a list depending on argument \code{type}.}
#'
#' \itemize{
#' \item type==makeMasterProblem: input is not used (empty list)
#' \item type==isProtectedSolution: input is a list of length 2 with elements 'input1' and 'input2'
#' \itemize{
#' \item element 'input1': numeric vector of calculated known lower cell bounds (from attacker's problem)
#' \item element 'input2': numeric vector of known upper cell bounds (from attacker's problem) }
#'
#' @return information from objects of class \code{problemInstance} depending on argument \code{type}
#' \itemize{
#' \item an object of class \code{linProb} if argument \code{type} matches 'makeMasterProblem'
#' \item logical vector of length 1 if argument \code{type} matches 'isProtectedSolution' with TRUE if all primary suppressed cells are adequately protected, FALSE otherwise }
#'
#' @export
#' @docType methods
#' @rdname calc.problemInstance-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.problemInstance', function(object, type, input) {standardGeneric('calc.problemInstance')})

# get-methods
setGeneric("g_sdcStatus", function(object) { standardGeneric("g_sdcStatus") })
setGeneric("g_primSupps", function(object) { standardGeneric("g_primSupps") })
setGeneric("g_secondSupps", function(object) { standardGeneric("g_secondSupps") })
setGeneric("g_forcedCells", function(object) { standardGeneric("g_forcedCells") })
setGeneric("g_type", function(object) { standardGeneric("g_type") })
setGeneric("g_freq", function(object) { standardGeneric("g_freq") })
setGeneric("g_strID", function(object) { standardGeneric("g_strID") })
setGeneric("g_UPL", function(object) { standardGeneric("g_UPL") })
setGeneric("g_LPL", function(object) { standardGeneric("g_LPL") })
setGeneric("g_SPL", function(object) { standardGeneric("g_SPL") })
setGeneric("g_nrVars", function(object) { standardGeneric("g_nrVars") })
setGeneric("g_lb", function(object) { standardGeneric("g_lb") })
setGeneric("g_ub", function(object) { standardGeneric("g_ub") })
setGeneric("g_w", function(object) { standardGeneric("g_w") })
setGeneric("g_numVars", function(object) { standardGeneric("g_numVars") })
setGeneric("g_hasPrimSupps", function(object) { standardGeneric("g_hasPrimSupps") })
setGeneric("g_hasSecondSupps", function(object) { standardGeneric("g_hasSecondSupps") })
setGeneric("g_hasForcedCells", function(object) { standardGeneric("g_hasForcedCells") })
setGeneric("g_weight", function(object) { standardGeneric("g_weight") })
setGeneric("g_suppPattern", function(object) { standardGeneric("g_suppPattern") })

# set-methods
setGeneric("s_sdcStatus<-", function(object, value) standardGeneric("s_sdcStatus<-"))
setGeneric("s_lb<-", function(object, value) standardGeneric("s_lb<-"))
setGeneric("s_ub<-", function(object, value) standardGeneric("s_ub<-"))
setGeneric("s_LPL<-", function(object, value) standardGeneric("s_LPL<-"))
setGeneric("s_UPL<-", function(object, value) standardGeneric("s_UPL<-"))
setGeneric("s_SPL<-", function(object, value) standardGeneric("s_SPL<-"))

# calc-methods
setGeneric("c_make_masterproblem", function(object, input) { standardGeneric("c_make_masterproblem") })
setGeneric("c_is_protected_solution", function(object, input) { standardGeneric("c_is_protected_solution") })
