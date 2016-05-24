## building_blocks.R
##   - Support for GP buildings blocks
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Support for GP buidling blocks
##'
##' Building blocks are a means for protecting expression subtrees
##' from modification through variation operators.
##' Often, certain functional units, represented as expression
##' subtrees in GP individuals, should stay intact during evolutionary
##' search. Building blocks at the leafs of expressions can be introduced
##' by adding them to the input variable set. Support for building blocks
##' is planned for a future release of RGP.
##'
##' \code{buildingBlock} transforms an R expression to a building block
##' to be used as an element of the input variable (or function) set. The
##' parameter \code{hardness} (a numerical value in the interval [0.0 , 1.0])
##' determines the protection strength against variation inside the building
##' building block. When \code{hardness} is set to \code{1.0} (the default),
##' the building block will never be subject to variaton through mutation
##' or crossover.
##' \code{buildingBlockq} is equivaltent to \code{buildingBlock}, but
##' quotes it's argument \code{expr} first.
##'
##' @param expr The expresion to transform to a building block.
##' @param hardness The strength of the protection against varition inside
##'   the building block. Must be a numeric in the interval [0.0, 1.0].
##'   A \code{hardness} of \code{1.0} (the default) means that the building
##'   block will never be subject to variation.
##' @return A building block.
##'
##' @rdname buildingBlocks 
##' @export
buildingBlock <- function(expr, hardness = 1.0)
  MapExpressionNodes(function(n) { buildingBlockTag(n) <- hardness; n },
                     expr, inners = TRUE, functions = TRUE, leafs = TRUE)

##' @rdname buildingBlocks 
##' @export
buildingBlockq <- function(expr, hardness = 1.0) buildingBlock(substitute(expr), hardness)

##' Building block tags
##'
##' To implement buidling blocks, i.e. subexpression protected from variation,
##' expression nodes may be tagged with \code{buildingBlockTags}. TODO
##' 
##' @param x An expression node.
##' @param value The value of the building block tag. Must be a numerical
##'   in the interval [0.0 1.0].
##'
##' @rdname buildingBlockTags
##' @export
buildingBlockTag <- function(x) {
  bbValue <- attr(x, "buildingBlockTag", exact = TRUE)
  if (is.null(bbValue)) 0.0 else bbValue
}

##' @rdname buildingBlockTags
##' @usage buildingBlockTag(x) <- value
##' @export
"buildingBlockTag<-" <- function(x, value) {
    attr(x, "buildingBlockTag") <- value
      x
}

##' @rdname buildingBlockTags
##' @export
hasBuildingBlockTag <- function(x) !is.null(buildingBlockTag(x))
