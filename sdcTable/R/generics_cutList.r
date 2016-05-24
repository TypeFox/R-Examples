#####################################
### Methods for class ''cutList'' ###
#####################################
#' query \code{cutList}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{cutList}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item constraints: constraint matrix of object
#' \item direction: directions of the constraints
#' \item rhs: right hand side of the constraints
#' \item nrConstraints: total number of constraints
#' @return information from objects of class \code{cutList} depending on argument \code{type}
#' \itemize{
#' \item object of class \code{simpleTriplet} if argument \code{type} matches 'constraints'
#' \item character vector if argument \code{type} matches 'direction'
#' \item numeric vector if argument \code{type} matches 'objective' or 'rhs'}
#'
#' @export
#' @docType methods
#' @rdname get.cutList-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("get.cutList", function(object, type) { standardGeneric("get.cutList")})

#' modify \code{cutList}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{cutList}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item addCompleteConstraint: add a constraint to argument \code{object}
#' \item removeCompleteConstraint: remove a constraint from argument \code{object}
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type==addCompleteConstraint: a list of length 1
#' \itemize{
#' \item first element: an object of class \code{cutList} with exactly one constraint }
#' \item type==removeCompleteConstraint: a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 specifying the index of the constraint that should be removed }
#'
#' @return an object of class \code{cutList}
#'
#' @export
#' @docType methods
#' @rdname set.cutList-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("set.cutList", function(object, type, input) { standardGeneric("set.cutList")})

#' perform calculations on \code{cutList}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{cutList}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item strengthen: strenghten constraints in argument \code{object}
#' \item checkViolation: check if a given solution violates any in argument \code{object}
#' \item bindTogether: combine two \code{cutList}-objects
#' @param input a list depending on argument \code{type}.}
#'
#' \itemize{
#' \item type==strengthen: input is not used (empty list)
#' \item type==checkViolation: input is a list of length 2
#' \itemize{
#' \item first element: numeric vector specifying a solution to a linear problem
#' \item second element: numeric vector specifying weights}
#' \item type==bindTogether: input is a list of length 1
#' \itemize{
#' \item first element: object of class \code{cutList} }
#'
#' @return manipulated data based on argument \code{type}
#' \itemize{
#' \item an object of class \code{cutList} if argument \code{type} matches 'strengthen' or 'bindTogether'
#' \item a logical vector of length 1 if argument \code{type} matches 'checkViolation' with TRUE if at least one constraint is violated by the given solution
#' }
#'
#' @export
#' @docType methods
#' @rdname calc.cutList-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("calc.cutList", function(object, type, input) { standardGeneric("calc.cutList")})

#' initialize \code{cutList}-objects depending on argument \code{type}
#'
#' @param type a character vector of length 1 defining what|how to initialize. Allowed types are:}
#' \itemize{
#' \item empty: create an empty \code{cutList}-object
#' \item singleCut: create a \code{cutList}-object with exactly one constraint
#' \item multipleCuts: create a \code{cutList}-object with more than one constraint
#' @param input a list depending on argument \code{type}.}
#'
#' \itemize{
#' \item type==empty: input is not used (empty list)
#' \item type==singleCut: input is a list of length 3
#' \itemize{
#' \item first element: numeric vector specifying a values for the row of the constraint matrix that must be created
#' \item second element: character vector of length 1 specifying the direction
#' \item third element: numeric vector of length 1 specifying the right hand side of the constraint}
#' \item type==multipleCuts: input is a list of length 3
#' \itemize{
#' \item first element: object of class \code{matrix}
#' \item second element: character vector specifying the direction of the constraints
#' \item third element: numeric vector specifying the right hand side of the constraints}
#'
#' @return an object of class \code{cutList}
#'
#' @export
#' @docType methods
#' @rdname init.cutList-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('init.cutList', function(type, input) { standardGeneric('init.cutList')})

# get-methods
setGeneric("g_constraints", function(object) { standardGeneric("g_constraints") })
setGeneric("g_direction", function(object) { standardGeneric("g_direction") })
setGeneric("g_rhs", function(object) { standardGeneric("g_rhs") })
setGeneric("g_nr_constraints", function(object) { standardGeneric("g_nr_constraints") })

# set-methods
setGeneric("s_add_complete_constraint<-", function(object, value) { standardGeneric("s_add_complete_constraint<-") })
setGeneric("s_remove_complete_constraint<-", function(object, value) { standardGeneric("s_remove_complete_constraint<-") })

# calc-methods
setGeneric("c_strengthen", function(object) { standardGeneric("c_strengthen") })
setGeneric("c_check_violation", function(object, input) { standardGeneric("c_check_violation") })
setGeneric("c_bind_together", function(object, input) { standardGeneric("c_bind_together") })


