#-----------------------------------#
# "PValueMat" CLASS DEFINITION    #
#-----------------------------------#
#' Univariate P-Value Matrix
#' 
#' Contains the raw \emph{p}-values and the \emph{p}-values adjusted for 
#' multiplicity for each pairwise comparison in each variable separately. 
#' This allows to see the contribution of each variable and each comparison 
#' to the final result
#' 
#' @title Class \code{PValueMat}
#' @aliases PValueMat
#' @aliases PValueMat-class
#' @aliases show,PValueMat,PValueMat-method
#' @aliases print,PValueMat,PValueMat-method
#' @rdname PValueMat-class
#' @keywords classes
#' @section Objects from the Class:
#'     Objects can be created by calls of the form 
#'     \code{new("PValueMat", ...)} but objects of the class are principally 
#'     supposed to be created and used internally for storing results of 
#'     \code{\link{NPC}}.
#' @section Slots:
#'     \describe{
#'     \item{\code{raw.p.values}:}{
#'         3-ways \code{array} containing the raw p.values}
#'     \item{\code{adj.p.values}:}{
#'         3-ways \code{array} containing the p.values adjusted for 
#'         multiplicity}
#'     \item{\code{p.adj.method}:}{
#'         multiplicity adjustment method employed}
#'     }
#' @section Prototype:
#'     \code{prototype} class has a \eqn{0 \times 0 \times 0}{0 x 0 x 0} 
#'     element for the first two slots and a \eqn{0}-length \code{character} 
#'     string
#' @section Methods:
#'     \describe{
#'     \item{show}{
#'         \code{signature(object = "PValueMat")}: shows only the main 
#'         information (on screen) for the object}
#'     \item{print}{
#'         \code{signature(x = "PValueMat")}: It prints the whole object 
#'         on screen (mostly useful for external saving)}
#'     }
#' @examples 
#'     showClass("PValueMat")
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @exportClass PValueMat

setClass("PValueMat",
    representation = representation(
        raw.p.values = "array",
        adj.p.values = "array",
        p.adj.method = "character"
    ),
    prototype = prototype(
        raw.p.values = array(NA, c(0, 0, 0)),
        adj.p.values = array(NA, c(0, 0, 0)),
        p.adj.method = character()
    )
)#END:Class

#---------------------------#
# "PValueMat" METHODS    #
#---------------------------#
##- show
#' \describe{
#' \item{\code{signature(object = "PValueMat")}}{
#'      Shows only the main information (on screen) for the object}
#' }
#' 
#' @name show
#' @aliases show,PValueMat-method
#' @docType methods
#' @rdname show-methods
#' @exportMethod show

setMethod(
    f = "show",
    signature = "PValueMat",
    definition = function(object) {
        cat("*** \"PValueMat\" object of package \"SOUP\" ***\n")
        ## raw.p.values p.values
        if(length(object@raw.p.values) > 0){
            cat("* Raw p.values *\n")
            print(round(object@raw.p.values, 5))
        } else { }
        ## adjusted p.values
        if(length(object@adj.p.values) > 0){
            cat("* Adjusted p.values, method: ", object@p.adj.method, " *\n", sep = "")
            print(round(object@adj.p.values, 5))
        } else { }
        cat("*** End \"PValueMat\" ***\n")
    }
)#END:Method

##- print
#' \describe{
#' \item{\code{signature(x = "PValueMat")}}{
#'      It prints the whole object on screen (mostly useful for 
#'      external saving)}
#' }
#' 
#' @name print
#' @docType methods
#' @aliases print,PValueMat-method
#' @rdname print-methods
#' @exportMethod print

setGeneric("print")
setMethod(
    f = "print", 
    signature = "PValueMat",
    definition = function(x, ...){
        cat("*** \"PValueMat\" object of package \"SOUP\" ***\n")
        ## raw.p.values p.values
        if(length(x@raw.p.values) > 0){
            cat("* Raw p.values *\n")
            print(x@raw.p.values)
        } else { }
        ## adjusted p.values
        if(length(x@adj.p.values) > 0){
            cat("* Adjusted p.values, method: ", x@p.adj.method, " *\n", sep = "")
            print(x@adj.p.values)
        } else { }
        cat("*** End \"PValueMat\" ***\n")
    }
)    

##- documenting it
# promptClass("PValueMat", file = "man/PValueMat.Rd")

