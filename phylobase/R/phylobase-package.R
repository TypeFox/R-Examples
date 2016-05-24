

##' Utilities and Tools for Phylogenetics
##'
##' Base package for phylogenetic structures and comparative data.
##'
##' \code{phylobase} provides a set of functions to associate and
##' manipulate phylogenetic information and data about the
##' species/individuals that are in the tree.
##'
##' \code{phylobase} intends to be robust, fast and efficient. We hope
##' other people use the data structure it provides to develop new
##' comparative methods in R.
##'
##' With \code{phylobase} it is easy to ensure that all your data are
##' represented and associated with the tips or the internal nodes of
##' your tree. \code{phylobase} provides functions to:
##' \itemize{
##'
##' \item prune (subset) your trees, find ancestor(s) a
##' descendant(s)
##'
##' \item find the most common recent ancestor of 2 nodes (MRCA)
##'
##' \item calculate the distance of a given node from the tip or
##' between two nodes in your tree
##'
##' \item robust functions to import data from NEXUS and Newick files
##' using the NEXUS Class Library (\url{https://github.com/mtholder/ncl/})
##' }
##'
##' @section History:
##'
##' \code{phylobase} was started during a Hackathlon at NESCent on
##' December 10-14 2007.
##'
##' Peter Cowan was a Google Summer of Code fellow in 2008 and
##' developed all the code for plotting.
##'
##' In December 2008, a mini-virtual Hackathlon was organized to clean
##' up and make the code more robust.
##'
##' In the spring and summer of 2009, Jim Regetz made several
##' contributions that made the code faster (in particular with the
##' re-ordering parts), found many bugs, and wrote most of the testing
##' code.
##'
##' \code{phylobase} was first released on CRAN on November 1st, 2009
##' with version 0.5.
##'
##' Since then, several releases have followed adding new
##' functionalities: better support of NEXUS files, creation of
##' \code{phylobase.options()} function that controls the \code{phylo4}
##' validator, rewrite of the validator in C++.
##'
##' Starting with 0.6.8, Francois Michonneau succeeds to Ben Bolker as
##' the maintainer of the package.
##'
##' @name phylobase-package
##' @aliases phylobase-package phylobase
##' @docType package
##' @section More Info:
##' See the help index \code{help(package="phylobase")} and run
##' \code{vignette("phylobase", "phylobase")} for further details and
##' examples about how to use \code{phylobase}.
##' @keywords package
##'
##' @useDynLib phylobase
##' @import methods
##' @import ape
##' @import RNeXML
##' @import grid
##' @import stats
##' @importFrom Rcpp evalCpp
##' @importFrom graphics plot
##' @importFrom utils head tail
##' @importFrom ade4 newick2phylog
##' @importFrom rncl rncl
##'
##' @exportMethod print head tail reorder plot summary
##'
## exportMethod should only be used for generics defined outside the package!
## @exportMethod phylo4 phylo4d
## @exportMethod edges edgeId hasEdgeLength edgeLength edgeLength<- sumEdgeLength edgeOrder
## @exportMethod isRooted rootNode rootNode<-
## @exportMethod isUltrametric
## @exportMethod subset prune [
## @exportMethod [<- [[ [[<-
## @exportMethod labels labels<- nodeLabels nodeLabels<- tipLabels tipLabels<- edgeLabels edgeLabels<-
## @exportMethod hasNodeLabels hasEdgeLabels hasDuplicatedLabels
NULL


##' Data from Darwin's finches
##'
##' Phylogenetic tree and morphological data for Darwin's finches, in different
##' formats
##'
##'
##' @name geospiza
##' @aliases geospiza geospiza_raw
##' @docType data
##' @format \code{geospiza} is a \code{phylo4d} object; \code{geospiza_raw} is a
##' list containing \code{tree}, a \code{phylo} object (the tree), \code{data},
##' and a data frame with the data (for showing examples of how to merge tree
##' and data)
##' @note Stolen from Luke Harmon's Geiger package, to avoid unnecessary
##' dependencies
##' @source Dolph Schluter via Luke Harmon
##' @keywords datasets
##' @examples
##'
##' data(geospiza)
##' plot(geospiza)
##'
NULL



##' 'Owls' data from ape
##'
##' A tiny tree, for testing/example purposes, using one of the examples from
##' the \code{ape} package
##'
##'
##' @name owls4
##' @docType data
##' @format This is the standard 'owls' tree from the \code{ape} package, in
##' \code{phylo4} format.
##' @source From various examples in the \code{ape} package
##' @keywords datasets
##' @examples
##'
##' data(owls4)
##' plot(owls4)
##'
NULL
