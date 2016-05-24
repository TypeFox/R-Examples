#' query \code{linProb}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{linProb}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item constraints: constraint matrix of object \code{linProb}
#' \item direction: directions of the constraints
#' \item rhs: right hand side of the constraints
#' \item objective: objective function
#' \item types: types of the objective variables
#' \item bounds: bounds of the objective variables
#' @return information from objects of class \code{linProb} depending on \code{type}
#' \itemize{
#' \item an object of class \code{simpleTriplet} if type matches 'constraints'
#' \item a character vector if type matches 'direction' or 'types'
#' \item a numeric vector if type matches 'objective' or 'rhs'
#' \item a list with elements 'lower' and 'upper' if type matches 'bounds'
#' \itemize{
#' \item element 'lower': a list with the first element containing indices and the second element containing corrsponding lower bounds
#' \item element 'upper': a list with the first element containing indices and the second element containing corrsponding upper bounds }
#' }
#'
#' @export
#' @docType methods
#' @rdname get.linProb-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('get.linProb', function(object, type) {standardGeneric('get.linProb')})

#' change \code{linProb}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{linProb}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item objective: change coefficients of the objective
#' \item direction: change vector of direction of the constraints
#' \item rhs: change vector of right hand side of the constraints
#' \item types: change vector of bounds of the objective variables
#' \item bounds: change bounds of the objective variables
#' \item constraints: change constraint matrix
#' \item removeCompleteConstraint: remove a specific constraint from the object
#' \item addCompleteConstraint: add a constraint to the object
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type==objective: a list of length 1
#' \itemize{
#' \item first element: numeric vector defining coefficients of the objective }
#' \item type==direction: a list of length 1
#' \itemize{
#' \item first element: character vector defining direction of the constraints }
#' \item type==rhs: a list of length 1
#' \itemize{
#' \item first element: numeric vector defining right hand side of the constraints }
#' \item type==types: a list of length 1
#' \itemize{
#' \item first element: character vector defining types of objective variables}
#' \item type==bounds: a list of length 2
#' \itemize{
#' \item element 'lower': a list with the first element containing indices and the second element containing corrsponding lower bounds
#' \item element 'upper': a list with the first element containing indices and the second element containing corrsponding upper bounds }
#' \item type==constraints: a list of length 1
#' \itemize{
#' \item first element: an object of class \code{simpleTriplet}}
#' \item type==removeCompleteConstraint: a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 defining the index of the constraint that should be removed }
#' \item type==addCompleteConstraint: a list of length 1
#' \itemize{
#' \item first element: an object of class \code{cutList} defining the constraint that should be added}
#' @return an object of class \code{linProb}
#'
#' @export
#' @docType methods
#' @rdname set.linProb-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('set.linProb', function(object, type, input) {standardGeneric('set.linProb')})

#' perform calculations on \code{linProb}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{linProb}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item solveProblem: solve the linear problem (minimize objective function)
#' \item fixVariables: try to fix objective variables to 0|1 based on dual costs depending on input
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type==solveProblem: a list of length 1
#' \itemize{
#' \item first element: character vector of length 1 specifying the solver to use. }
#' \item type==fixVariables: a list of length 3
#' \itemize{
#' \item first element: numeric vector specifying lower bounds for the objective variables
#' \item second element: numeric vector specifying upper bounds for the objective variables
#' \item third element: numeric vector specifying indices of primary suppressed cells }
#' @return manipulated data based on argument \code{type}
#' \itemize{
#' \item list containing the solution and additional information if argument \code{type} matches 'solveProblem
#' \item a numeric vector of indices if argument \code{type} matches 'fixVariables'
#' }
#'
#' @export
#' @docType methods
#' @rdname calc.linProb-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.linProb', function(object, type, input) {standardGeneric('calc.linProb')})

# get-methods
setGeneric("g_constraints", function(object) { standardGeneric("g_constraints") })
setGeneric("g_direction", function(object) { standardGeneric("g_direction") })
setGeneric("g_rhs", function(object) { standardGeneric("g_rhs") })
setGeneric("g_objective", function(object) { standardGeneric("g_objective") })
setGeneric("g_types", function(object) { standardGeneric("g_types") })
setGeneric("g_bounds", function(object) { standardGeneric("g_bounds") })

# set-methods
setGeneric("s_objective<-", function(object, value) { standardGeneric("s_objective<-") })
setGeneric("s_direction<-", function(object, value) { standardGeneric("s_direction<-") })
setGeneric("s_rhs<-", function(object, value) { standardGeneric("s_rhs<-") })
setGeneric("s_types<-", function(object, value) { standardGeneric("s_types<-") })
setGeneric("s_bounds<-", function(object, value) { standardGeneric("s_bounds<-") })
setGeneric("s_constraints<-", function(object, value) { standardGeneric("s_constraints<-") })
setGeneric("s_add_complete_constraint<-", function(object, value) { standardGeneric("s_add_complete_constraint<-") })
setGeneric("s_remove_complete_constraint<-", function(object, value) { standardGeneric("s_remove_complete_constraint<-") })

# calc-methods
setGeneric("c_solve_problem", function(object, input) { standardGeneric("c_solve_problem") })
setGeneric("c_fix_variables", function(object, input) { standardGeneric("c_fix_variables") })
