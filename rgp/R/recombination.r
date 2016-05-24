## recombination.R
##   - Functions for recombining GP individuals (e.g. by crossover)
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Select random childs or subtrees of an expression
##'
##' \code{randchild} returns a uniformly random direct child of an expression.
##' \code{randsubtree} returns a uniformly random subtree of an expression. Note that this subtree
##' must not be a direct child.
##
##' @param expr The expression to select random childs or subtrees from.
##' @param subtreeprob The probability for \code{randsubtree} to select a certain subtree instead
##' of searching further via an recursive call.
##'
##' @rdname randomExpressionChilds
##' @export
randchild <- function(expr)
  expr[[ceiling(runif(1, 1, length(expr)))]]

##' @rdname randomExpressionChilds
##' @export
randsubtree <- function(expr, subtreeprob = 0.1)
  if (is.call(expr) && runif(1) > subtreeprob) randsubtree(randchild(expr), subtreeprob) else expr

##' Random crossover (recombination) of functions and expressions
##'
##' Replace a random subtree of \code{func1} (\code{expr1}) with a random subtree of
##' \code{func2} (\code{expr2}) and return the resulting function (expression), i.e.
##' the modified \code{func1} (\code{expr1}).
##' \code{crossoverexpr} handles crossover of expressions instead of functions.
##' \code{crossoverexprFast} is a fast (i.e. implemented in efficient C code)
##' albeit less flexible variant of \code{crossoverexpr}.
##' \code{crossoverTyped} and \code{crossoverexprTyped} only exchage replace subtress
##' if the sTypes of their root nodes match.
##' \code{crossoverTwoPoint} is a variant of \code{crossover} that swaps subtrees
##' at uniform randomly selected points and returns both children.
##' \code{crossoverexprTwoPoint} works analogously for expressions.
##'
##' All RGP recombination operators operating on functions have the S3 class
##' \code{c("recombinationOperator", "function")}.
##'
##' @param expr1 The first parent R expression.
##' @param func1 The first parent R function.
##' @param expr2 The second parent R expression.
##' @param func2 The second parent R function.
##' @param crossoverprob The probability of crossover at each node of the first parent
##'   function (expression).
##' @param breedingFitness A breeding function. See the documentation for
##'   \code{\link{geneticProgramming}} for details.
##' @param breedingTries The number of breeding steps.
##' @return The child function (expression) or functions (expressions).
##'
##' @rdname expressionCrossover
##' @export
crossover <- function(func1, func2, crossoverprob = 0.1,
                      breedingFitness = function(individual) TRUE,
                      breedingTries = 50) {
  doCrossover <- function() {
    child <- new.function(envir = environment(func1)) 
    formals(child) <- formals(func1)
    body(child) <- crossoverexpr(body(func1), body(func2), crossoverprob)
    child
  }
  breed(doCrossover, breedingFitness, breedingTries)
}
class(crossover) <- c("recombinationOperator", "function")

##' @rdname expressionCrossover
##' @export
crossoverexpr <- function(expr1, expr2, crossoverprob) {
  newexpr1 <- expr1
  if (is.call(expr1)) {
    crossoverindexA <- ceiling(runif(1, 1, length(expr1)))
    if (runif(1) <= crossoverprob) { # try to do crossover here
      if (runif(1) > buildingBlockTag(newexpr1[[crossoverindexA]]))
        newexpr1[[crossoverindexA]] <- randsubtree(expr2, crossoverprob)
    } else { # try to do crossover somewhere below in this tree
      newexpr1[[crossoverindexA]] <- crossoverexpr(expr1[[crossoverindexA]], expr2, crossoverprob)
    }
  }
  newexpr1
}

##' @rdname expressionCrossover
##' @export
crossoverexprFast <- function(expr1, expr2)
  .Call("crossover_single_point_R", expr1, expr2)[[1]]

##' @rdname expressionCrossover
##' @export
crossoverexprTwoPoint <- function(expr1, expr2) {
  index1 <- floor(runif(1, 0, exprSize(expr1)))
  index2 <- floor(runif(1, 0, exprSize(expr2)))
  child1 <- replaceSubtreeAt(expr1, index1, subtreeAt(expr2, index2))
  child2 <- replaceSubtreeAt(expr2, index2, subtreeAt(expr1, index1))
  list(child1, child2)  
}

##' @rdname expressionCrossover
##' @export
crossoverTyped <- function(func1, func2, crossoverprob = 0.1,
                           breedingFitness = function(individual) TRUE,
                           breedingTries = 50) {
  doCrossoverTyped <- function() {
    child <- new.function(envir = environment(func1)) 
    formals(child) <- formals(func1)
    body(child) <- crossoverexprTyped(body(func1), body(func2), crossoverprob)
    child
  }
  breed(doCrossoverTyped, breedingFitness, breedingTries)
}
class(crossoverTyped) <- c("recombinationOperator", "function")

##' @rdname expressionCrossover
##' @export
crossoverexprTyped <- function(expr1, expr2, crossoverprob) {
  newexpr1 <- expr1
  if (is.call(expr1)) {
    crossoverindex <- ceiling(runif(1, 1, length(expr1)))
    if (runif(1) <= crossoverprob) { # try to do crossover here
      if (runif(1) > buildingBlockTag(newexpr1[[crossoverindex]])) {
        ## select a random subtree of expr2 as a replacement expression...
        replacementExpression <- randsubtree(expr2, crossoverprob)
        ## collect the range types of the crossover points...
        crossoverExpressionRangeType <- rangeTypeOfType(sType(newexpr1[[crossoverindex]]))
        replacementExpressionRangeType <- rangeTypeOfType(sType(replacementExpression))
        ## only do crossover if the range types of the crossover points match...
        if (identical(replacementExpressionRangeType, crossoverExpressionRangeType)) {
          newexpr1[[crossoverindex]] <- replacementExpression
        }
      }
    } else { # try to do crossover somewhere below in this tree
      newexpr1[[crossoverindex]] <- crossoverexprTyped(expr1[[crossoverindex]], expr2, crossoverprob)
    }
  }
  newexpr1
}
