#' perform calculations on multiple objects depending on argument \code{type}
#'
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item makePartitions: information on subtables required for HITAS and HYPECUBE algorithms
#' \item genMatMFull: the constraint matrix used in the master problem
#' \item makeAttackerProblem: set up the attackers problem for a given (sub)table
#' \item calcFullProblem: calculate a complete problem object containing all information required to solve the secondary cell suppression problem
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item if type matches 'makePartitions', 'genMatMFull' or 'makeAttackerProblem': a list of length 2 with elements 'objectA' and 'objectB'
#' \itemize{
#' \item element 'object A': an object of class \code{problemInstance}
#' \item element 'object B': an object of class \code{dimInfo} }
#' \item type matches 'calcFullProblem': a list of length 1
#' \itemize{
#' \item element 'object A': an object of class \code{dataObj}
#' \item element 'object B': an object of class \code{dimInfo} }
#'
#' @return manipulated data based on argument \code{type}
#' \itemize{
#' \item list with elements 'groups', 'indices', 'strIDs', 'nrGroups' and 'nrTables' if argument \code{type} matches 'makePartitions'
#' \item object of class \code{simpleTriplet} if argument \code{type} matches 'genMatMFull'
#' \item object of class \code{linProb} if argument \code{type} matches 'makeAttackerProblem'
#' \item object of class \code{sdcProblem} if argument \code{type} matches 'calcFullProblem'
#' }
#'
#' @export
#' @docType methods
#' @rdname calc.multiple-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.multiple', function(type, input) {standardGeneric('calc.multiple')})

setGeneric("c_make_partitions", function(input) { standardGeneric("c_make_partitions") })
setGeneric("c_gen_mat_m", function(input) { standardGeneric("c_gen_mat_m") })
setGeneric("c_make_att_prob", function(input) { standardGeneric("c_make_att_prob") })
setGeneric("c_calc_full_prob", function(input) { standardGeneric("c_calc_full_prob") })
