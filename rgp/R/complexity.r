## complexity.R
##   - Complexity measures for GP individuals and estimators of GP search space size
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Complexity measures for R functions and expressions
##'
##' \code{exprDepth} returns the depth of the tree representation ("exression tree") of an R expression.
##' \code{funcDepth} returns the tree depth of the body expression of an R function.
##' \code{exprSize} returns the number of nodes in the tree of an R expression.
##' \code{exprLeaves} returns the number of leave nodes in the tree of an R expression.
##' \code{exprCount} returns the number of tree nodes in an R expression matching a given predicate.
##' \code{funcSize} returns the number of nodes in the body expression tree of an R function.
##' \code{funcLeaves} returns the number of leave nodes in the body expression tree of an R function.
##' \code{funcCount} returns the number of nodes in an R function body expression matching a given predicate.
##' \code{exprVisitationLength} returns the visitation length of the tree of an R expression.
##' The visitation length is the total number of nodes in all possible subtrees of a tree.
##' \code{funcVisitationLength} returns the visitation length of the body expression tree of an R function.
##' \code{fastExprVisitationLength} and \code{fastFuncVisitationLength} are variants written in optimized
##' C code.
##' The visitation length can be interpreted as the size of the expression obtained by substituting all
##' inner functions by their function bodies (see "Crossover Bias in Genetic Programming", Maarten Keijzer
##' and James Foster).
##'
##' @param expr An R expression.
##' @param func An R function.
##' @param predicate An R predicate (function with range type \code{logical}).
##' @param intermediateResults Whether to return complexity measures for all subtrees also.
##'
##' @rdname expressionComplexityMeasures
##' @export
exprDepth <- function(expr)
  if (is.call(expr)) {
    max(as.vector(Map(exprDepth, rest(as.list(expr))), mode = "integer")) + 1
  } else 1

##' @rdname expressionComplexityMeasures
##' @export
funcDepth <- function(func) exprDepth(body(func))

##' @rdname expressionComplexityMeasures
##' @export
exprSize <- function(expr)
  if (is.call(expr)) {
    sum(as.vector(Map(exprSize, rest(as.list(expr))), mode = "integer")) + 1
  } else 1

##' @rdname expressionComplexityMeasures
##' @export
exprLeaves <- function(expr)
  if (is.call(expr)) {
    sum(as.vector(Map(exprLeaves, rest(as.list(expr))), mode = "integer"))
  } else 1

##' @rdname expressionComplexityMeasures
##' @export
exprCount <- function(expr, predicate = function(node) TRUE) {
  isCall <- is.call(expr)
  matchesPredicate <- predicate(expr)
  if (isCall && matchesPredicate) {
    sum(as.vector(Map(function(e) exprCount(e, predicate = predicate), rest(as.list(expr))), mode = "integer")) + 1
  } else if (isCall) {
    sum(as.vector(Map(function(e) exprCount(e, predicate = predicate), rest(as.list(expr))), mode = "integer"))
  } else if (matchesPredicate) 1 else 0
}

##' @rdname expressionComplexityMeasures
##' @export
funcSize <- function(func) exprSize(body(func))

##' @rdname expressionComplexityMeasures
##' @export
funcLeaves <- function(func) exprLeaves(body(func))

##' @rdname expressionComplexityMeasures
##' @export
funcCount <- function(func, predicate = function(node) TRUE) exprCount(body(func), predicate = predicate)

##' @rdname expressionComplexityMeasures
##' @export
exprVisitationLength <- function(expr, intermediateResults = FALSE) {
  results <- c()
  exprVisitationLengthRecursive <- function(expr) {
    ## The visitation length of a tree T is the sum of the number of nodes of all subtrees of T...
    if (is.call(expr)) {
      childrenResults <- lapply(rest(expr), exprVisitationLengthRecursive)
      childrenSums <- Reduce(function(a,b) c(a[1] + b[1], a[2] + b[2]), childrenResults, c(0,0))
      childrenSizesSum <- childrenSums[1]
      childrenVisitationLengthsSum <- childrenSums[2]
      thisTreeSize <- 1 + childrenSizesSum
      thisTreeVisitationLength <- thisTreeSize + childrenVisitationLengthsSum
      if (intermediateResults) results <<- c(thisTreeVisitationLength, results)
      c(thisTreeSize, thisTreeVisitationLength)
    } else {
      if (intermediateResults) results <<- c(1, results)
      c(1, 1)
    }
  }
  if (intermediateResults) {
    exprVisitationLengthRecursive(expr)
    results
  } else exprVisitationLengthRecursive(expr)[[2]]
}

##' @rdname expressionComplexityMeasures
##' @export
fastExprVisitationLength <- function(expr, intermediateResults = FALSE)
  .Call("sexp_visitation_length_R", expr, intermediateResults) 

##' @rdname expressionComplexityMeasures
##' @export
funcVisitationLength <- function(func, intermediateResults = FALSE)
  exprVisitationLength(body(func), intermediateResults = intermediateResults)

##' @rdname expressionComplexityMeasures
##' @export
fastFuncVisitationLength <- function(func, intermediateResults = FALSE)
  .Call("func_visitation_length_R", func, intermediateResults) 

##' Upper bounds for expression tree search space sizes
##'
##' These functions return the number of structurally different expressions or expression shapes of a given
##' depth or size that can be build from a fixed function- and input-variable set. Here, "expression shape"
##' means the shape of an expression tree, not taking any node labels into account.
##' \code{exprShapesOfDepth} returns the number of structurally different expression shapes of a depth
##' exactly equal to \code{n}.
##' \code{exprShapesOfMaxDepth} returns the number of structurally different expression shapes of a depth
##' less or equal to \code{n}.
##' \code{exprsOfDepth} returns the number of structurally different expressions of a depth exactly equal to
##' \code{n}. Note that constants are handled by conceptually substiting them with a fresh input variable.
##' \code{exprShapesOfMaxDepth} returns the number of structurally different expressions of a depth
##' less or equal to \code{n}. Note that constants are handled by conceptually substiting them with a fresh
##' input variable.
##' \code{exprShapesOfSize}, \code{exprShapesOfMaxSize}, \code{exprsOfSize}, \code{exprsOfMaxSize} are
##' equivalents that regard expression tree size (number of nodes) instead of expression tree depth.
##'
##' @param funcset The function set.
##' @param inset The set of input variables.
##' @param n The fixed size or depth.
##'
##' @rdname expressionCounts
##' @export
exprShapesOfDepth <- function(funcset, n)
  if (n < 1) {
    0
  } else if (n == 1) {
    1 # there is only one expression shape of depth 1
  } else {
    exprShapesOfDepthNminusOne <- exprShapesOfDepth(funcset, n - 1)
    uniqueArities <- unique(funcset$arities)
    exprsWithFixedRootArity <- exprShapesOfDepthNminusOne ^ uniqueArities
    sum(exprsWithFixedRootArity)
  }

##' @rdname expressionCounts
##' @export
exprShapesOfMaxDepth <- function(funcset, n)
  if (n < 1) {
    0
  } else if (n == 1) {
    1 # there is only one expression shape of max depth 1
  } else {
    exprShapesOfDepthNminusOne <- exprShapesOfMaxDepth(funcset, n - 1)
    uniqueArities <- unique(funcset$arities)
    exprsWithFixedRootArity <- exprShapesOfDepthNminusOne ^ uniqueArities
    sum(exprsWithFixedRootArity) + 1
  }

##' @rdname expressionCounts
##' @export
exprsOfDepth <- function(funcset, inset, n)
  if (n < 1) {
    0
  } else if (n == 1) {
    1 + length(inset) # an expression of depth 1 is either a constant or an input variable
  } else {
    exprsOfDepthNminusOne <- exprsOfDepth(funcset, inset, n - 1)
    arities <- funcset$arities
    exprsWithFixedRoot <- exprsOfDepthNminusOne ^ arities
    sum(exprsWithFixedRoot)
  }

##' @rdname expressionCounts
##' @export
exprsOfMaxDepth <- function(funcset, inset, n)
  if (n < 1) {
    0
  } else if (n == 1) {
    1 + length(inset) # an expression of depth 1 is either a constant or an input variable
  } else {
    exprsOfMaxDepthNminusOne <- exprsOfMaxDepth(funcset, inset, n - 1)
    arities <- funcset$arities
    exprsWithFixedRoot <- exprsOfMaxDepthNminusOne ^ arities
    sum(exprsWithFixedRoot) + 1 + length(inset)
  }

##' @rdname expressionCounts
##' @export
exprShapesOfSize <- function(funcset, n) NA # TODO

##' @rdname expressionCounts
##' @export
exprShapesOfMaxSize <- function(funcset, n) NA # TODO

##' @rdname expressionCounts
##' @export
exprsOfSize <- function(funcset, inset, n) NA # TODO

##' @rdname expressionCounts
##' @export
exprsOfMaxSize <- function(funcset, inset, n) NA # TODO
