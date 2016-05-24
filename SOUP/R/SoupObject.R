#-------------------------------------#
#    "SoupObject" CLASS DEFINITION    #
#-------------------------------------#
setClassUnion("callOrNull", c("NULL", "call"))
setClassUnion("listOrRankResults", c("RankResults", "arrayOrList"))

#' SOUP Main Object
#' 
#' Contains the main results of the \pkg{SOUP} analysis
#' 
#' @title Class \code{SoupObject}
#' @aliases SoupObject
#' @aliases SoupObject-class
#' @aliases show,SoupObject,SoupObject-method
#' @aliases print,SoupObject,SoupObject-method
#' @aliases plot,SoupObject,SoupObject-method
#' @rdname SoupObject-class
#' @keywords classes
#' @section Objects from the Class:
#'     Objects can be created by calls of the form 
#'     \code{new("SoupObject", ...)}. 
#'     It contains all results from the analysis performed by the 
#'     \code{\link{SOUP}} function.
#' @section Slots:
#'     \describe{
#'     \item{\code{call}:}{
#'         a \code{call} object that contains the call of the 
#'         \code{\link{SOUP}} main function}
#'     \item{\code{pValueMat}:}{
#'         \code{\linkS4class{PValueMat}} object containing univariate 
#'         \emph{p}-values}
#'     \item{\code{rankResults}:}{
#'         \code{\linkS4class{RankResults}} object containing the results of 
#'         the ranking}
#'     \item{\code{permSpace}:}{
#'         \code{\linkS4class{PermSpace}} object containing the permutation 
#'         spaces}
#'     }
#' @section Methods:
#'     \describe{
#'     \item{show}{
#'         \code{signature(object = "SoupObject")}: shows only the main 
#'         information (on screen) for the object}
#'     \item{print}{
#'         \code{signature(x = "SoupObject")}: It prints the whole object 
#'         on screen (mostly useful for external saving)}
#'     }
#' @examples 
#'     showClass("SoupObject")
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @exportClass SoupObject


##- SoupObject
setClass("SoupObject",
    representation = representation(
        call         = "callOrNull",
        pValueMat   = "PValueMat",
        rankResults  = "listOrRankResults",
        permSpace    = "PermSpace"
    ),
    prototype = prototype(
        call         = NULL,
        pValueMat   = new("PValueMat"),
        rankResults  = new("RankResults"),
        permSpace    = new("PermSpace")
    )    
)#-END-

#---------------------------#
#    'SoupObject' METHODS    #
#---------------------------#
##- show
#' Methods for function \code{show} in package \pkg{SOUP}
#' 
#' \describe{
#' \item{\code{signature(object = "SoupObject")}}{
#'      Shows only the main information (on screen) for the object}
#' }
#' 
#' @name show
#' @aliases show,SoupObject-method
#' @docType methods
#' @rdname show-methods
#' @exportMethod show

setMethod(
    f = "show",
    signature = "SoupObject",
    definition = function(object) {
        call4print <- as.character(enquote(object@call))[2L]
        call4print <- strsplit(call4print, split = ",")[[1]]
        
        cat("*** \"SoupObject\" object of package \"SOUP\" ***\n")
        cat("* Call:\n") 
        cat(call4print, sep = ", ", fill = TRUE)
        cat("\n\n")
        # show(object@pValueMat)
        if(is.list(object@rankResults))
        {
            SubsWtsNames <- names(object@permSpace@T.H0Low)
            for(indRun in seq_along(object@rankResults))
            {
                cat("\n")
                cat("Subsets And Weights: ", SubsWtsNames[indRun], 
                        fill = TRUE)
                cat("\n")
                show(object@rankResults[[indRun]])
                cat("\n")
            }# END: for
        } else {
            show(object@rankResults)
        }# END:ifelse-list
        # if(is.list(object@permSpace)) {
            # lapply(object@permSpace, FUN = show)
        # } else {        
            # show(object@permSpace)
        # }# END:ifelse-list
        cat("* Seed for the RNG: ", object@permSpace@seed, "*\n")
        cat("*** End of \"SoupObject\" ***\n")
    }
)


##- print
#' Methods for function \code{show} in package \pkg{SOUP}
#' 
#' \describe{
#' \item{\code{signature(x = "SoupObject")}}{
#'      It prints the whole object on screen (mostly useful for 
#'      external saving)}
#' }
#' 
#' @name print
#' @docType methods
#' @aliases print,SoupObject-method
#' @rdname print-methods
#' @exportMethod print

setGeneric("print")
setMethod(
    f = 'print', 
    signature = "SoupObject",
    definition = function(x, ...) {
        cat("*** \"SoupObject\" object of package \"SOUP\" ***\n")
        cat("* Call:\n")
        cat(as.character(enquote(x@call))[2L], "\n")
        print(x@pValueMat)
        if(is.list(x@rankResults)) {
            lapply(x@rankResults, FUN = print)
        } else {
            print(x@rankResults)
        }# END:ifelse-list
        if(is.list(x@permSpace)) {
            lapply(x@permSpace, FUN = print)
        } else {        
            print(x@permSpace)
        }# END:ifelse-list
        cat("* Seed for the RNG: ", x@permSpace@seed, "*\n")
        cat("*** End of \"SoupObject\" ***\n")
    }
)

##- documenting it
# promptClass("SoupObject", file = "man/SoupObjecte.Rd")
