##' The phylo4 class
##' 
##' Classes for phylogenetic trees
##' 
##' @name phylo4-class
##' @docType class
##' @section Objects from the Class: Phylogenetic tree objects can be created by
##' calls to the \code{\link{phylo4}} constructor function.  Translation
##' functions from other phylogenetic packages are also available. See
##' \code{\link{coerce-methods}}.
##' @author Ben Bolker, Thibaut Jombart
##' @seealso The \code{\link{phylo4-methods}} constructor, the
##' \code{\link{checkPhylo4}} function to check the validity of
##' \code{phylo4} objects. See also the \code{\link{phylo4d-methods}}
##' constructor and the \linkS4class{phylo4d} class.
##' @keywords classes
##' @include RcppExports.R checkdata.R
##' @export
setClass("phylo4",
         representation(edge = "matrix",
                        edge.length = "numeric",
                        label = "character",
                        edge.label = "character",
                        order = "character",
                        annote = "list"),
         prototype = list(
                        edge = matrix(nrow = 0, ncol = 2,
                            dimname = list(NULL, c("ancestor", "descendant"))),
                        edge.length = numeric(0),
                        label = character(0),
                        edge.label = character(0),
                        order = "unknown",
                        annote = list()
                       ),
         validity = checkPhylo4)
