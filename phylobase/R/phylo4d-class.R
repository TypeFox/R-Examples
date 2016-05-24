###################################
## phylo4d class
## extend: phylo with data
##' phylo4d class
##' 
##' S4 class for phylogenetic tree and data.
##' 
##' 
##' @name phylo4d-class
##' @docType class
##' @section Objects from the Class: Objects can be created from various trees
##' and a data.frame using the constructor \code{phylo4d}, or using
##' \code{new("phylo4d", \dots{})} for empty objects.
##' @author Ben Bolker, Thibaut Jombart
##' @seealso \code{\link{coerce-methods}} for translation
##' functions. The \code{\link{phylo4d-methods}} constructor. See also
##' the \code{\link{phylo4-methods}} constructor, the
##' \linkS4class{phylo4} class, and the \code{\link{checkPhylo4}}
##' function to check the validity of \code{phylo4} trees.
##' @keywords classes
##' @export
##' @include phylo4-methods.R formatData.R
##' @examples
##'   example(read.tree, "ape")
##'   obj <- phylo4d(as(tree.owls.bis,"phylo4"), data.frame(wing=1:3))
##'   obj
##'   names(obj)
##'   summary(obj)
setClass("phylo4d",
         representation(data="data.frame",
                        metadata = "list"),

         prototype = list(
           data = data.frame(NULL),
           metadata = list()),

         validity = checkPhylo4,
         contains = "phylo4")

