#' query \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item dataObj: a list containing the (raw) input data
#' \item problemInstance: return the current problem instance
#' \item partition: a list containing information on the subtables that are required to be protected as well as information on the processing order of the subtables
#' \item elapsedTime: the elapsed time of the protection algorithm so far
#' \item dimInfo: information on the variables defining the hierarchical table
#' \item indicesDealtWith: a set of indices that have already been dealt with during the protection algorithmus
#' \item startI: current level at which subtables need to be protected (useful when restarting HITAS|HYPERCUBE)
#' \item startJ: current number of the subtable within a given level that needs to be protected (useful when restarting HITAS|HYPERCUBE)
#' \item innerAndMarginalCellInfo: for a given problem, get indices of inner- and marginal table cells
#'
#' @return information from objects of class \code{sdcProblem} depending on argument \code{type}
#' \itemize{
#' \item an object of class \code{dataObj} (or NULL) if \code{type} matches 'dataObj'
#' \item an object of class \code{problemInstance} (or NULL) if \code{type} matches 'problemInstance'
#' \item a list (or NULL) if argument \code{type} matches 'partition' containing the following elements:
#' \itemize{
#' \item element 'groups': list with each list-element being a character vector specifying a specific level-group
#' \item element 'indices': list with each list-element being a numeric vector defining indices of a subtable
#' \item element 'strIDs': list with each list-element being a character vector defining IDs of a subtable
#' \item element 'nrGroups': numeric vector of length 1 defining the total number of groups that have to be considered
#' \item element 'nrTables': numeric vector of length 1 defining the total number of subtables that have to be considered}
#' \item a list (or NULL) if argument \code{type} matches 'innerAndMarginalCellInfo' containing the following elements:
#' \itemize{
#' \item element 'innerCells': character vector specifying ID's of inner cells
#' \item element 'totCells': character vector specifying ID's of marginal cells
#' \item element 'indexInnerCells': numeric vector specifying indices of inner cells
#' \item element 'indexTotCells': numeric vector specifying indices of marginal cells}
#' \item an object of class \code{dimInfo} (or NULL) if \code{type} matches 'dimInfo'
#' \item numeric vector if argument \code{type} matches 'elapsedTime'
#' \item numeric vector of length 1 if argument \code{type} matches 'startI' or 'startJ'
#' }
#'
#' @export
#' @docType methods
#' @rdname get.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('get.sdcProblem', function(object, type) {standardGeneric('get.sdcProblem')})

#' modify \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item problemInstance: set|modify slot 'problemInstance' of argument \code{object}
#' \item partition: set|modify slot 'partition' of argument \code{object}
#' \item startI: set|modify slot 'startI' of argument \code{object}
#' \item startJ: set|modify slot 'startJ' of argument \code{object}
#' \item indicesDealtWith: set|modify slot 'indicesDealtWith' of argument \code{object}
#' \item elapsedTime: set|modify slot 'elapsedTime' of argument \code{object}
#' @param input a list with elements depending on argument \code{type}.}
#'
#' \itemize{
#' \item an object of class \code{problemInstance} if argument \code{type} matches 'problemInstance'
#' \item a list (derived from calc.multiple(type='makePartition', ...) if argument \code{type} matches 'partition'
#' \item a numeric vector of length 1 if argument \code{type} matches 'startI', 'startJ' or 'elapsedTime'
#' \item a numeric vector if argument \code{type} matches 'indicesDealtWith'
#'
#' @return an object of class \code{sdcProblem}
#'
#' @export
#' @docType methods
#' @rdname set.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('set.sdcProblem', function(object, type, input) {standardGeneric('set.sdcProblem')})

#' perform calculations on \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item rule.freq: modify suppression status within \code{object} according to frequency suppression rule
#' \item rule.nk: modify sdcStatus of \code{object} according to nk-dominance rule
#' \item rule.p: modify sdcStatus of \code{object} according to p-percent rule
#' \item rule.pq: modify sdcStatus of \code{object} according to pq-rule
#' \item heuristicSolution: obtain a heuristic (greedy) solution to the problem defined by \code{object}
#' \item cutAndBranch: solve a secondary cell suppression problem defined by \code{object} using cut and branch
#' \item anonWorker: is used to solve the suppression problem depending on information provided with argument \code{input}
#' \item ghmiter: solve a secondary cell suppression problem defined by \code{object} using hypercube algorithm
#' \item preprocess: perform a preprocess procedure by trying to identify primary suppressed cells that are already protected due to other primary suppressed cells
#' \item cellID: find index of cell defined by information provided with argument \code{input}
#' \item finalize: create an object of class \code{safeObj}
#' \item ghmiter.diagObj: calculate codes required to identify diagonal cells given a valid cell code - used for ghmiter-algorithm only
#' \item ghmiter.calcInformation: calculate information for quaders identified by diagonal indices - used for ghmiter-algorithm only
#' \item ghmiter.suppressQuader: suppress a quader based on indices
#' \item ghmiter.selectQuader: select a quader for suppression depending on information provided with argument \code{input} - used for ghmiter-algorithm only
#' \item ghmiter.suppressAdditionalQuader: select and suppress an additional quader (if required) based on information provided with argument \code{input} - used for ghmiter-algorithm only
#' \item contributingIndices: calculate indices within the current problem that contribute to a given cell
#' \item reduceProblem: reduce the problem given by \code{object} using a vector of indices
#' \item genStructuralCuts: calculate cuts that are absolute necessary for a valid solution of the secondary cell suppression problem
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item a list (typically generated using genParaObj()) specifying parameters for primary cell suppression if argument \code{type} matches 'rule.freq', 'rule.nk' or 'rule.p'
#' \item a list if argument \code{type} matches 'heuristicSolution' having the following elements:
#' \itemize{
#' \item element 'aProb': an object of class \code{linProb} defining the attacker's problem
#' \item element 'validCuts': an object of class \code{cutList} representing a list of constraints
#' \item element 'solver': a character vector of length 1 specifying a solver to use
#' \item element 'verbose': a logical vector of length 1 setting if verbose output is desired }
#' \item a list (typically generated using genParaObj()) specifying parameters for the secondary cell suppression problem if argument \code{type} matches 'cutAndBranch', 'anonWorker', 'ghmiter', 'preprocess'
#' \item a list of length 3 if argument \code{type} matches 'cellID' having following elements
#' \itemize{
#' \item first element: character vector specifying variable names that need to exist in slot 'dimInfo' of \code{object}
#' \item second element: character vector specifying codes for each variable that define a specific table cell
#' \item third element: logical vector of length 1 with TRUE setting verbosity and FALSE to turn verbose output off}
#' \item a list of length 3 if argument \code{type} matches 'ghmiter.diagObj' having following elements
#' \itemize{
#' \item first element: numeric vector of length 1
#' \item second element:  a list with as many elements as dimensional variables have been specified and each element being a character vector of dimension-variable specific codes
#' \item third element: logical vector of length 1 defining if diagonal indices with frequency == 0 should be allowed or not}
#' \item a list of length 4 if argument \code{type} matches 'ghmiter.calcInformation' having following elements
#' \itemize{
#' \item first element: a list object typically generated with method \code{calc.sdcProblem} and type=='ghmiter.diagObj'
#' \item second element: a list with as many elements as dimensional variables have been specified and each element being a character vector of dimension-variable specific codes
#' \item third element: numeric vector of length 1 specifying a desired protection level
#' \item fourth element: logical vector of length 1 defining if quader containing empty cells should be allowed or not}
#' \item a list of length 1 if argument \code{type} matches 'ghmiter.suppressQuader' having following element
#' \itemize{
#' \item first element: numeric vector of indices that should be suppressed }
#' \item a list of length 2 if argument \code{type} matches 'ghmiter.selectQuader' having following elements
#' \itemize{
#' \item first element: a list object typically generated with method \code{calc.sdcProblem} and type=='ghmiter.calcInformation'
#' \item second element: a list (typically generated using genParaObj())}
#' \item a list of length 4 if argument \code{type} matches 'ghmiter.suppressAdditionalQuader' having following elements
#' \itemize{
#' \item first element: a list object typically generated with method \code{calc.sdcProblem} and type=='ghmiter.diagObj'
#' \item second element: a list object typically generated with method \code{calc.sdcProblem} and type=='ghmiter.calcInformation'
#' \item third element: a list object typically generated with method \code{calc.sdcProblem} and type=='ghmiter.selectQuader'
#' \item fourth element: a list (typically generated using genParaObj()) }
#' \item a list of length 1 if argument \code{type} matches 'contributingIndices' having following element
#' \itemize{
#' \item first element: character vector of length 1 being an ID for which contributing indices should be calculated }
#' \item a list of length 1 if argument \code{type} matches 'reduceProblem' having following element
#' \itemize{
#' \item first element: numeric vector defining indices of cells that should be kept in the reduced problem }
#' \item an empty list if argument \code{type} matches 'genStructuralCuts'
#' @return information from objects of class \code{sdcProblem} depending on argument \code{type}
#' \itemize{
#' \item an object of class \code{sdcProblem} if argument \code{type} matches 'rule.freq', 'rule.nk', 'rule.p', 'cutAndBranch', 'anonWorker', 'ghmiter', 'ghmiter.supressQuader', 'ghmiter.suppressAdditionalQuader' or 'reduceProblem'
#' \item a numeric vector with elements being 0 or 1 if argument \code{type} matches 'heuristicSolution'
#' \item a list if argument \code{type} matches 'preprocess' having following elements:
#' \itemize{
#' \item element 'sdcProblem': an object of class \code{sdcProblem}
#' \item element 'aProb': an object of class \code{linProb}
#' \item element 'validCuts': an object of class \code{cutList}	}
#' \item a numeric vector of length 1 specifying the index of the cell of interest if argument \code{type} matches 'cellID'
#' \item an object of class \code{safeObj} if argument \code{type} matches 'finalize'
#' \item a list if argument \code{type} matches 'ghmiter.diagObj' having following elements:
#' \itemize{
#' \item element 'cellToProtect': character vector of length 1 defining the ID of the cell to protect
#' \item element 'indToProtect': numeric vector of length 1 defining the index of the cell to protect
#' \item element 'diagIndices': numeric vector defining indices of possible cells defining cubes }
#' \item a list containing information about each quader that could possibly be suppressed if argument \code{type} matches 'ghmiter.calcInformation'
#' \item a list containing information about a single quader that should be suppressed if argument \code{type} matches 'ghmiter.selectQuader'
#' \item a numeric vector with indices that contribute to the desired table cell if argument \code{type} matches 'contributingIndices'
#' \item an object of class \code{cutList} if argument \code{type} matches 'genStructuralCuts'
#' }
#' @export
#' @docType methods
#' @rdname calc.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.sdcProblem', function(object, type, input) {standardGeneric('calc.sdcProblem')})

# get-methods
setGeneric("g_problemInstance", function(object) { standardGeneric("g_problemInstance") })
setGeneric("g_dimInfo", function(object) { standardGeneric("g_dimInfo") })
setGeneric("g_partition", function(object) { standardGeneric("g_partition") })
setGeneric("g_elapsedTime", function(object) { standardGeneric("g_elapsedTime") })
setGeneric("g_dataObj", function(object) { standardGeneric("g_dataObj") })
setGeneric("g_startI", function(object) { standardGeneric("g_startI") })
setGeneric("g_startJ", function(object) { standardGeneric("g_startJ") })
setGeneric("g_indicesDealtWith", function(object) { standardGeneric("g_indicesDealtWith") })
setGeneric("g_innerAndMarginalCellInfo", function(object) { standardGeneric("g_innerAndMarginalCellInfo") })
setGeneric("g_df", function(object, ...) { standardGeneric("g_df") })

# set-methods
setGeneric("s_problemInstance<-", function(object, value) standardGeneric("s_problemInstance<-"))
setGeneric("s_partition<-", function(object, value) standardGeneric("s_partition<-"))
setGeneric("s_startI<-", function(object, value) standardGeneric("s_startI<-"))
setGeneric("s_startJ<-", function(object, value) standardGeneric("s_startJ<-"))
setGeneric("s_indicesDealtWith<-", function(object, value) standardGeneric("s_indicesDealtWith<-"))
setGeneric("s_elapsedTime<-", function(object, value) standardGeneric("s_elapsedTime<-"))

# calc-methods
setGeneric("c_rule_freq", function(object, input) { standardGeneric("c_rule_freq") })
setGeneric("c_rule_nk", function(object, input) { standardGeneric("c_rule_nk") })
setGeneric("c_rule_nk", function(object, input) { standardGeneric("c_rule_nk") })
setGeneric("c_rule_p", function(object, input) { standardGeneric("c_rule_p") })
setGeneric("c_rule_pq", function(object, input) { standardGeneric("c_rule_pq") })
setGeneric("c_heuristic_solution", function(object, input) { standardGeneric("c_heuristic_solution") })
setGeneric("c_anon_worker", function(object, input) { standardGeneric("c_anon_worker") })
setGeneric("c_opt_cpp", function(object, input) { standardGeneric("c_opt_cpp") })
setGeneric("c_hitas_cpp", function(object, input) { standardGeneric("c_hitas_cpp") })
setGeneric("c_quick_suppression", function(object, input) { standardGeneric("c_quick_suppression") })
setGeneric("c_cut_and_branch", function(object, input) { standardGeneric("c_cut_and_branch") })
setGeneric("c_ghmiter", function(object, input) { standardGeneric("c_ghmiter") })
setGeneric("c_preprocess", function(object, input) { standardGeneric("c_preprocess") })
setGeneric("c_cellID", function(object, input) { standardGeneric("c_cellID") })
setGeneric("c_finalize", function(object, input) { standardGeneric("c_finalize") })
setGeneric("c_ghmiter_diag_obj", function(object, input) { standardGeneric("c_ghmiter_diag_obj") })
setGeneric("c_ghmiter_calc_info", function(object, input) { standardGeneric("c_ghmiter_calc_info") })
setGeneric("c_ghmiter_suppress_quader", function(object, input) { standardGeneric("c_ghmiter_suppress_quader") })
setGeneric("c_ghmiter_select_quader", function(object, input) { standardGeneric("c_ghmiter_select_quader") })
setGeneric("c_ghmiter_supp_additional", function(object, input) { standardGeneric("c_ghmiter_supp_additional") })
setGeneric("c_contributing_indices", function(object, input) { standardGeneric("c_contributing_indices") })
setGeneric("c_reduce_problem", function(object, input) { standardGeneric("c_reduce_problem") })
setGeneric("c_gen_structcuts", function(object, input) { standardGeneric("c_gen_structcuts") })


