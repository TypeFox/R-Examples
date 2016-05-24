#-----------------------------------#
# "RankResults" CLASS DEFINITION    #
#-----------------------------------#
#' Contains results of the algorithm, if \code{n} values of \code{alpha} are 
#' provided for to the \pkg{SOUP} function then there will be \code{n} 
#' (possibly different) rankings, one for each value of \code{alpha}
#' 
#' @title Ranking Results
#' @aliases RankResults RankResults-class
#' @rdname RankResults-class
#' @keywords classes
#' @section Objects from the Class:
#'     Objects can be created by calls of the form 
#'     \code{new("RankResults", ...)}. 
#' @section Slots:
#'     \describe{
#'     \item{\code{alpha}:}{
#'         \code{numeric} vector containing the values of \code{alpha} used 
#'         for rejecting each pairwise hypothesis}
#'     \item{\code{ranks}:}{
#'         \code{matrix} containing the rankings obtained from the ranking 
#'         algorithm}
#'     \item{\code{p.values}:}{
#'         \code{matrix} containing the p-values used for the rankings}
#'     \item{\code{p.adj.method}:}{
#'         \code{character} string indicating which multiplicity adjustment 
#'         method was used for the pairwise p-values}
#'     }
#' @section Prototype:
#'     \code{prototype} class has a \eqn{0 \times 0 \times 0}{0 x 0 x 0} 
#'     element for the first two slots and a \eqn{0}-length \code{character} 
#'     string
#' @section Methods:
#'     \describe{
#'     \item{show}{
#'         \code{signature(object = "RankResults")}: shows only the main 
#'         information (on screen) for the object}
#'     \item{print}{
#'         \code{signature(x = "RankResults")}: It prints the whole object 
#'         on screen (mostly useful for external saving)}
#'     \item{as.list}{
#'      Coerce the \code{RankResults} object to a list of \code{RankResults}
#'      objects}
#'     }
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @examples 
#'     showClass("RankResults")
#' @exportClass RankResults
#' 
setClass("RankResults",
    representation = representation(
        alpha        = "numeric",
        ranks        = "matrix",
        p.values     = "matrix",
        p.adj.method = "character"
    ),
    prototype = prototype(
        alpha        = numeric(),
        ranks        = matrix(NA, 0, 0),
        p.values     = matrix(NA, 0, 0),
        p.adj.method = character()
    )
)#END:Class



#---------------------------#
# "RankResults" METHODS    #
#---------------------------#
#' Methods for function \code{as.list} in package \pkg{SOUP}
#' 
#' @name as.list
#' @aliases as.list,RankResults-method
#' @docType methods
#' @rdname RankResults-method
#' @section Methods:
#'     \describe{
#'     \item{\code{signature(object = "RankResults")}}{
#'          conversion of the object (\code{RankResults}) into a \code{list}}
#'     }
#' @exportMethod as.list
setGeneric("as.list", function(x) standardGeneric("as.list"))

setMethod(
    f = "as.list", 
    signature = "RankResults",
    definition = function(x, ...)
    {
        if (length(x) > 1L)
        {
            as.list.numeric_version(x, ...)
        } else
        {
            return(x)
        }
    }
)#END:Method



##- show
##' Method show
##' @name show
##' @rdname show-methods
##' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @name show
#' @aliases show,RankResults-method
#' @docType methods
#' @rdname show-methods
#' @exportMethod show
#setGeneric("show", function(x) standardGeneric("show"))
setMethod(
    f = "show", 
    signature = "RankResults",
    definition = function(object) {
        cat("*** \"RankResults\" object of package \"SOUP\" ***\n")
        K <- ncol(object@p.values) * (ncol(object@p.values) - 1) / 2
        if(length(object@p.values) > 1) {
            cat("* Final ranking results *\n")
            # cat(paste("* alpha =", paste(object@alpha, collapse = ", "), "*", collapse = ""), "\n")
            print(object@ranks)
            cat("\n* Associated p.value matrix *\n")
            print(round(object@p.values, 5))
            cat("\n* p.values multiplicity adjustment method:", 
                    object@p.adj.method, "*\n")
            # cat("*** End \"RankResults\" ***\n")
        } else { }
    }
)#END:Method

##- print
#' @name print
#' @docType methods
#' @aliases print,RankResults-method
#' @rdname print-methods
#' @exportMethod print
setGeneric("print", function(x) standardGeneric("print"))
setMethod(
    f = "print", 
    signature = "RankResults",
    definition = function(x, ...) {
        if(length(x@p.values) > 1){
            cat("*** \"RankResults\" object of package \"SOUP\" ***\n")
            K <- ncol(x@p.values) * (ncol(x@p.values) - 1) / 2
            cat("* Final ranking results *\n")
            cat(paste("* alpha =", paste(x@alpha, collapse = ", "), "*", collapse = ""), "\n")
            print(x@ranks)
            cat("* Associated p.value matrix *\n")
            print(x@p.values)
            cat("* p.values multiplicity adjustment method:", x@p.adj.method, "*\n")
            # cat("*** End \"RankResults\" ***\n")
        } else { }
    }
)#END:Method

