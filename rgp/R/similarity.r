## similarity.R
##   - Similarity and distance measures for R functions and expressions
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Patrick Koch, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Similarity and Distance Measures for R Functions and Expressions
##'
##' These functions implement several similarity and distance measures for R functions
##' (i.e. their body expressions).
##' TODO check and document measure-theoretic properties of each measure defined here
##' TODO these distance measures are metrics, some of them are norm-induced metrics
##' \code{commonSubexpressions} returns the set of common subexpressions of \code{expr1}
##' and \code{expr2}. This is not a metric by itself, but can be used to implement
##' several subtree-based similarity metrics.
## '\code{numberOfcommonSubexpressions} returns the number of common subexpressions
##' of \code{expr1} and \code{expr2}.
##' \code{sizeWeightedNumberOfcommonSubexpressions} returns the number of common
##' subexpressions of \code{expr1} and \code{expr2}, weighting the size of each common
##' subexpression. Note that for every expression \emph{e},
##' \code{sizeWeightedNumberOfcommonSubexpressions(} \emph{e} \code{, } \emph{e}
##' \code{) == exprVisitationLength(} \emph{e} \code{)}.
##' \code{normalizedNumberOfCommonSubexpressions} returns the ratio of the number of
##' common subexpressions of \code{expr1} and \code{expr2} in relation to the number
##' of subexpression in the larger expression of \code{expr1} and \code{expr2}.
##' \code{normalizedSizeWeightedNumberOfcommonSubexpressions} returns the ratio of
##' the size-weighted number of common subexpressions of \code{expr1} and \code{expr2}
##' in relation to the visitation length of the larger expression of \code{expr1} and
##' \code{expr2}.
##' \code{NCSdist} and \code{SNCSdist} are distance metrics derived from
##' \code{normalizedNumberOfCommonSubexpressions} and
##' \code{normalizedSizeWeightedNumberOfCommonSubexpressions} respectively.
##' \code{differingSubexpressions}, and code{numberOfDifferingSubexpressions}
##' are duals of the functions described above, based on counting the number of
##' differing subexpressions of \code{expr1} and \code{expr2}. The possible functions
##' "normalizedNumberOfDifferingSubexpressions" and
##' "normalizedSizeWeightedNumberOfDifferingSubexpressions" where ommited because they
##' are always equal to \code{NCSdist} and \code{SNCSdist} by definition.
##' \code{trivialMetric} The "trivial" metric M(a, b) that is 0 iff a == b, 1 otherwise.
##' \code{normInducedTreeDistance} Uses a norm on expression trees and a metric on tree
##' node labels to induce a metric M on expression trees A and B: If both A and B are empty
##' (represented as \code{NULL}), M(A, B) := 0. If exactly one of A or B is empty, M(A, B) :=
##' "the norm applied to the non-empty tree". If neither A or B is empty, the difference
##' of their root node labels (as measured by \code{labelDistance}) is added to the sum of
##' the differences of the children. The children lists are padded with empty trees to
##' equalize their sizes. The summation operator can be changed via \code{distanceFoldOperator}.
##' \code{normInducedFunctionDistance} Is wrapper that applies \code{normInducedTreeDistance}
##' to the bodies of the given functions.
##'
##' @param expr1 An R expression.
##' @param expr2 An R expression.
##' @param a An R object.
##' @param b An R object.
##' @param norm A norm to derive a tree distance metric from.
##' @param labelDistance A metric for measuring distances of tree node labels, i.e. function
##'   names or constants.
##' @param distanceFoldOperator The operator used by \code{normInducedTreeDistance} to combine
##'   the measures subtree distances, defaults to `+`.
##'
##' @rdname expressionSimilarityMeasures
##' @export
commonSubexpressions <- function(expr1, expr2)
  intersect(subexpressions(expr1), subexpressions(expr2))

##' @rdname expressionSimilarityMeasures
##' @export
numberOfCommonSubexpressions <- function(expr1, expr2)
  length(commonSubexpressions(expr1, expr2))

##' @rdname expressionSimilarityMeasures
##' @export
normalizedNumberOfCommonSubexpressions <- function(expr1, expr2)
  numberOfCommonSubexpressions(expr1, expr2) / max(exprSize(expr1), exprSize(expr2))

##' @rdname expressionSimilarityMeasures
##' @export
NCSdist <- function(expr1, expr2)
  1 - normalizedNumberOfCommonSubexpressions(expr1, expr2)

##' @rdname expressionSimilarityMeasures
##' @export
sizeWeightedNumberOfCommonSubexpressions <- function(expr1, expr2) {
  weightedCommonSubexpressions <- sapply(commonSubexpressions(expr1, expr2), exprSize)
  if (identical(weightedCommonSubexpressions, list())) 0 else sum(weightedCommonSubexpressions)
}

##' @rdname expressionSimilarityMeasures
##' @export
normalizedSizeWeightedNumberOfCommonSubexpressions <- function(expr1, expr2) {
  largerExpr <- if (exprSize(expr1) >= exprSize(expr2)) expr1 else expr2
  sizeWeightedNumberOfCommonSubexpressions(expr1, expr2) / exprVisitationLength(largerExpr)
}

##' @rdname expressionSimilarityMeasures
##' @export
SNCSdist <- function(expr1, expr2)
  1 - normalizedSizeWeightedNumberOfCommonSubexpressions(expr1, expr2)

##' @rdname expressionSimilarityMeasures
##' @export
differingSubexpressions <- function(expr1, expr2) {
  expr1subs <- subexpressions(expr1)
  expr2subs <- subexpressions(expr2)
  setdiff(union(expr1subs, expr2subs), intersect(expr1subs, expr2subs))
}

##' @rdname expressionSimilarityMeasures
##' @export
numberOfDifferingSubexpressions <- function(expr1, expr2)
  length(differingSubexpressions(expr1, expr2))

##' @rdname expressionSimilarityMeasures
##' @export
sizeWeightedNumberOfDifferingSubexpressions <- function(expr1, expr2) {
  weightedDifferingSubexpressions <- sapply(differingSubexpressions(expr1, expr2), exprSize)
  if (identical(weightedDifferingSubexpressions, list())) 0 else sum(weightedDifferingSubexpressions)
}

##' @rdname expressionSimilarityMeasures
##' @export
trivialMetric <- function(a, b) if (a == b) 0 else 1

##' @rdname expressionSimilarityMeasures
##' @export
normInducedTreeDistance <- function(norm, labelDistance = trivialMetric, distanceFoldOperator = NULL) {
  if (is.null(distanceFoldOperator)) distanceFoldOperator <- `+` # ...to make Roxygen happy
  M <- function(expr1, expr2)
    if (is.null(expr1) && is.null(expr2)) { # both trees are empty
      0
    } else if (is.null(expr1)) { # one tree is empty, apply the norm
      norm(expr2)
    } else if (is.null(expr2)) { # - " -
      norm(expr1)
    } else { # none of the trees are empty, recurse
      labelDistance <- trivialMetric(exprLabel(expr1), exprLabel(expr2))
      childrenOfExpr1 <- exprChildrenOrEmptyList(expr1)
      childrenOfExpr2 <- exprChildrenOrEmptyList(expr2)
      maxExprWidth <- max(length(childrenOfExpr1), length(childrenOfExpr2))
      # "pad" the children lists with "NULL" to equalize their lengths
      nullPaddedChildrenOfExpr1 <-
        append(childrenOfExpr1, replicate(maxExprWidth - length(childrenOfExpr1), NULL))
      nullPaddedChildrenOfExpr2 <-
        append(childrenOfExpr2, replicate(maxExprWidth - length(childrenOfExpr2), NULL))
      childrenDistances <- Map(M, nullPaddedChildrenOfExpr1, nullPaddedChildrenOfExpr2)
      distanceFoldOperator(labelDistance,
                           Reduce(distanceFoldOperator, childrenDistances, 0))
    }
  M
}

##' @rdname expressionSimilarityMeasures
##' @export
normInducedFunctionDistance <- function(norm, labelDistance = trivialMetric, distanceFoldOperator = NULL) {
  treeM <- normInducedTreeDistance(norm, labelDistance, distanceFoldOperator)
  function(fun1, fun2) treeM(body(fun1), body(fun2))
}

##' A \code{dist} function that supports custom metrics
##'
##' This function computes and returns the distance matrix computed by using the given metric
##' to compute the distances between the rows of a data list or vector. Note that in contrast
##' to \code{\link{dist}}, \code{x} has to be a vector and the the distance \code{metric}
##' is an arbitrary function that must be symmetric and definite.
##' @seealso \code{\link{dist}}
##'
##' @param x A vector or list of objects.
##' @param metric A metric, i.e. a function of two arguments that returns a numeric. Note
##' that a metric must be definite and symmetric, otherwise the results will be undefined.
##' @param diag \code{TRUE} iff the diagonal of the distance matrix should be printed by
##' \code{print.dist}.
##' @param upper \code{TRUE} iff the upper triangle of the distance matrix should be printed
##' by \code{print.dist}.
##' @return A distance matrix.
##'
##' @rdname customDist
##' @export
customDist <- function(x, metric, diag = FALSE, upper = FALSE) {
  xlen <- length(x)
  dm <- matrix(nrow = xlen, ncol = xlen)
  for (j in 1:(xlen-1)) for (i in j:(xlen-1)) dm[i+1,j] <- metric(x[[i+1]], x[[j]])
  as.dist(dm, upper = upper, diag = diag)
}

##' Return the Children of an Expression or the Empty List if there are None
##'
##' Internal tool function that returns the children expressions of an R expression or the
##' empty list if there are no children, i.e. if the expression is atomic or \code{NULL}.
##' If the expression is a "function" expression, i.e. an expression that would evaluate
##' to a function, \code{exprChildrenOrEmptyList} will return the function body expression
##' as the only child.
##'
##' @param expr The expression to return the children for.
##' @return The expression's children as a list, or the empty list if there are none.
exprChildrenOrEmptyList <- function(expr)
  if (inherits(expr, "call") && expr[[1]] == as.name("function")) { # a function expression
    list(expr[[3]]) # third element is the "body" of the function expression
  } else if (inherits(expr, "call")) { # a call
    as.list(expr[-1])
  } else list() # an atom

##' Return the "label" at the Root Node of an Expression Tree
##'
##' Internal tool function that returns the function name if \code{expr} is a call,
##' or otherwise just \code{expr} itself.
##'
##' @param expr The expression to return the root label for.
##' @return The expression's root label.
exprLabel <- function(expr) if (inherits(expr, "call")) expr[[1]] else expr
