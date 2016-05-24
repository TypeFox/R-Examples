## mutation.R
##   - Functions for mutating GP individuals
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Random mutation of functions and expressions
##'
##' RGP implements two sets of mutation operators. The first set is inspired by classical
##' GP systems. Mutation strength is controlled by giving mutation probabilities:
##' \code{mutateFunc} mutates a function \eqn{f} by recursively replacing inner function labels in
##'   \eqn{f} with probability \code{mutatefuncprob}.
##' \code{mutateSubtree} mutates a function by recursively replacing inner nodes with
##'   newly grown subtrees of maximum depth \code{maxsubtreedepth}.
##' \code{mutateNumericConst} mutates a function by perturbing each numeric (double) constant \eqn{c}
##'   with probability \code{mutateconstprob} by setting \eqn{c := c + rnorm(1, mean = mu, sd = sigma)}.
##'   Note that constants of other typed than \code{double} (e.g \code{integer}s) are not affected.
##'
##' \code{mutateFuncTyped}, \code{mutateSubtreeTyped}, and \code{mutateNumericConstTyped} are
##' variants of the above functions that only create well-typed result expressions.
##'
##' \code{mutateFuncFast}, \code{mutateSubtreeFast}, \code{mutateNumericConstFast} are variants
##' of the above untyped mutation function implemented in C. They offer a considerably faster
##' execution speed for the price of limited flexibility. These variants take function bodies
##' as arguments (obtain these via R's \code{body} function) and return function bodies as results.
##' To turn a function body into a function, use RGP's \code{\link{makeClosure}} tool function.
##'
##' The second set of mutation operators features a more orthogonal design, with each individual
##' operator having a only a small effect on the genotype. Mutation strength is controlled by
##' the integral \code{strength} parameter.
##' \code{mutateChangeLabel} Selects a node (inner node or leaf) by uniform random sampling and replaces
##'   the label of this node by a new label of matching type.
##' \code{mutateInsertSubtree} Selects a leaf by uniform random sampling and replaces it with a matching
##'   subtree of the exact depth of \code{subtreeDepth}.
##' \code{mutateDeleteSubtree} Selects a subree of the exact depth of \code{subtreeDepth} by uniform random
##'   sampling and replaces it with a matching leaf.
##' \code{mutateChangeDeleteInsert} Either applies \code{mutateChangeLabel}, \code{mutateInsertSubtree},
##'   or \code{mutateDeleteSubtree}. The probability weights for selecting an operator can be supplied
##'   via the ...Probability arguments (probability weights are normalized to a sum of 1). 
##' \code{mutateDeleteInsert} Either applies \code{mutateDeleteSubtree} or \code{mutateInsertSubtree}. The
##'  probability weights for selecting an operator can be supplied via the ...Probability arguments
##'  (probability weights are normalized to a sum of 1).
##' The above functions automatically create well-typed result expressions when used in a strongly
##' typed GP run.
##'
##' All RGP mutation operators have the S3 class \code{c("mutationOperator", "function")}.
##'
##' @param func The function to mutate randomly.
##' @param funcbody The function body to mutate randomly, obtain it via \code{body(func)}.
##' @param funcset The function set.
##' @param inset The set of input variables.
##' @param conset The set of constant factories.
##' @param mutatefuncprob The probability of trying to replace an inner function at each node.
##' @param mutatesubtreeprob The probability of replacing a subtree with a newly grown subtree
##'   at each node.
##' @param maxsubtreedepth The maximum depth of newly grown subtrees.
##' @param mutateconstprob The probability of mutating a constant by adding \code{rnorm(1)} to it.
##' @param strength The number of individual point mutations (changes, insertions, deletions) to
##'   perform.
##' @param subtreeDepth The depth of the subtrees to insert or delete.
##' @param constprob The probability of creating a constant versus an input variable.
##' @param insertprob The probability to insert a subtree.
##' @param deleteprob The probability to insert a subtree.
##' @param constmin The lower limit for numeric constants.
##' @param constmax The upper limit for numeric onstants.
##' @param mu The normal distribution mean for random numeric constant mutation.
##' @param sigma The normal distribution standard deviation for random numeric constant mutation.
##' @param subtreeprob The probability of creating a subtree instead of a leaf in the random subtree
##'   generator function.
##' @param iterations The number of times to apply a mutation operator to a GP individual. This
##'   can be used as a generic way of controling the strength of the genotypic effect of mutation. 
##' @param changeProbability The probability for selecting the \code{mutateChangeLabel} operator.
##' @param deleteProbability The probability for selecting the \code{mutateDeleteSubtree} operator.
##' @param insertProbability The probability for selecting the \code{mutateInsertSubtree} operator.
##' @param breedingFitness A breeding function. See the documentation for
##'   \code{\link{geneticProgramming}} for details.
##' @param breedingTries The number of breeding steps.
##' @return The randomly mutated function.
##'
##' @rdname expressionMutation
##' @export
mutateFunc <- function(func, funcset, mutatefuncprob = 0.1,
                       breedingFitness = function(individual) TRUE,
                       breedingTries = 50) {
  mutatefuncexpr <- function(expr, funcset, mutatefuncprob) {
    if (is.call(expr)) {
      oldfunc <- expr[[1]]
      if (runif(1) > buildingBlockTag(oldfunc)) {
        newfunccandidate <- if (runif(1) <= mutatefuncprob)
            toName(randelt(funcset$all, prob = attr(funcset$all, "probabilityWeight")))
          else oldfunc
        newfunc <- if (arity(newfunccandidate) == arity(oldfunc)) newfunccandidate else oldfunc
        as.call(append(newfunc, Map(function(e) mutatefuncexpr(e, funcset, mutatefuncprob), rest(expr))))
      } else expr
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutatefuncexpr(body(func), funcset, mutatefuncprob)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateFunc) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateSubtree <- function(func, funcset, inset, conset, mutatesubtreeprob = 0.1, maxsubtreedepth = 5,
                          breedingFitness = function(individual) TRUE,
                          breedingTries = 50) {
  mutatesubtreeexpr <- function(expr, funcset, inset, conset, mutatesubtreeprob, maxsubtreedepth) {
    if (runif(1) <= mutatesubtreeprob) { # replace current node with new random subtree
      if (runif(1) > buildingBlockTag(expr)) {
        randexprGrow(funcset, inset, conset, maxdepth = maxsubtreedepth)
      } else expr
    } else if (is.call(expr)) {
      as.call(append(expr[[1]],
                     Map(function(e) mutatesubtreeexpr(e, funcset, inset, conset,
                                                       mutatesubtreeprob, maxsubtreedepth),
                         rest(expr))))
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutatesubtreeexpr(body(func), funcset, inset, conset, mutatesubtreeprob, maxsubtreedepth)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateSubtree) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateNumericConst <- function(func, mutateconstprob = 0.1,
                               breedingFitness = function(individual) TRUE,
                               breedingTries = 50, mu = 0.0, sigma = 1.0) {
  mutateconstexpr <- function(expr, mutateconstprob) {
    if (is.call(expr)) {
      as.call(append(expr[[1]], Map(function(e) mutateconstexpr(e, mutateconstprob), rest(expr))))
    } else if (runif(1) <= mutateconstprob && is.double(expr)) {
      if (runif(1) > buildingBlockTag(expr)) {
        mutatedExpr <- expr + rnorm(1, mean = mu, sd = sigma)
        mutatedExpr
      } else expr
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = environment(func))
    formals(mutant) <- formals(func)
    body(mutant) <- mutateconstexpr(body(func), mutateconstprob)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateNumericConst) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateFuncTyped <- function(func, funcset, mutatefuncprob = 0.1,
                            breedingFitness = function(individual) TRUE,
                            breedingTries = 50) {
  mutatefuncexprTyped <- function(expr, funcset, mutatefuncprob) {
    if (is.call(expr)) { # only look at calls, this mutation ignores terminal nodes...
      oldfunc <- expr[[1]]
      if (runif(1) > buildingBlockTag(oldfunc)) {
        oldfuncType <- sType(oldfunc)
        oldfuncRangeType <- rangeTypeOfType(oldfuncType)
        ## Select a candidate for a new function of matching range type. This can of course result
        ## in a candidate function with a different domain type. If this happens the mutation is
        ## simply aborted, because searching again for a matching function would cost too much time...
        newfunccandidate <- if (runif(1) <= mutatefuncprob)
          toName(randelt(funcset$byRange[[oldfuncRangeType$string]],
                         prob = attr(funcset$byRange[[oldfuncRangeType$string]], "probabilityWeight")))
        else oldfunc
        newfunccandidateType <- sType(newfunccandidate)
        newfunc <- if (identical(newfunccandidateType, oldfuncType)) newfunccandidate else oldfunc
        newcall <- as.call(append(newfunc, Map(function(e) mutatefuncexprTyped(e, funcset, mutatefuncprob), rest(expr))))
        newcall
      } else expr
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutatefuncexprTyped(body(func), funcset, mutatefuncprob)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateFuncTyped) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateSubtreeTyped <- function(func, funcset, inset, conset, mutatesubtreeprob = 0.1, maxsubtreedepth = 5,
                               breedingFitness = function(individual) TRUE,
                               breedingTries = 50) {
  mutatesubtreeexprTyped <- function(expr, funcset, inset, conset, mutatesubtreeprob, maxsubtreedepth) {
    if (runif(1) <= mutatesubtreeprob) { # replace current node with new random subtree of correct type
      if (runif(1) > buildingBlockTag(expr)) {
        type <- rangeTypeOfType(sType(expr))
        randexprTypedGrow(type, funcset, inset, conset, maxdepth = maxsubtreedepth)
      } else expr
    } else if (is.call(expr)) {
      mutatedExpr <-
        as.call(append(expr[[1]],
                       Map(function(e) mutatesubtreeexprTyped(e, funcset, inset, conset,
                                                              mutatesubtreeprob, maxsubtreedepth),
                           rest(expr))))
      mutatedExpr
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutatesubtreeexprTyped(body(func), funcset, inset, conset, mutatesubtreeprob, maxsubtreedepth)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateSubtreeTyped) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateNumericConstTyped <- function(func, mutateconstprob = 0.1,
                                    breedingFitness = function(individual) TRUE,
                                    breedingTries = 50) {
  mutateconstexprTyped <- function(expr, mutateconstprob) {
    if (is.call(expr)) {
      mutatedExpr <- as.call(append(expr[[1]], Map(function(e) mutateconstexprTyped(e, mutateconstprob), rest(expr))))
      mutatedExpr
    } else if (runif(1) <= mutateconstprob && is.double(expr)) {
      if (runif(1) > buildingBlockTag(expr)) {
        mutatedExpr <- expr + rnorm(1)
        mutatedExpr
      } else expr
    } else expr
  }
  doMutation <- function() {
    mutant <- new.function(envir = environment(func))
    formals(mutant) <- formals(func)
    body(mutant) <- mutateconstexprTyped(body(func), mutateconstprob)
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateNumericConstTyped) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateChangeLabel <- function(func, funcset, inset, conset,
                              strength = 1,
                              breedingFitness = function(individual) TRUE,
                              breedingTries = 50) {
  numberOfNodes <- funcSize(func)
  sampledMutationPoints <- sample.int(numberOfNodes, replace = FALSE,
                                      size = min(strength, numberOfNodes))
  currentNode <- 0 # here, local mutable state is more efficient than recursion
  mutateExpressionChangeLabel <- function(expr) {
    currentNode <<- currentNode + 1
    if (is.symbol(expr)) {
      if (currentNode %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
        if (hasStype(expr)) {
          newInputVariable <- toName(randelt(inset$byType[[sType(expr)$string]],
                                             prob = attr(inset$byType[[sType(expr)$string]], "probabilityWeight")))
          newInputVariable
        } else {
          newInputVariable <- toName(randelt(inset$all, prob = attr(inset$all, "probabilityWeight")))
          newInputVariable
        }
      } else expr
    } else if (is.call(expr)) {
      if (identical(expr[[1]], as.symbol("function"))) {
        stop("mutateChnageLabel: Support for anonymous function nodes is not implemented.") # TODO
      } else if (identical(expr[[1]], as.symbol("("))) {
        ## Just skip parentheses in the expression tree...
        restExpr <- rest(expr)
        mutatedExpr <- as.call(append(expr[[1]], Map(mutateExpressionChangeLabel, restExpr)))
        mutatedExpr
      } else {
        mutatedLabel <- if (currentNode %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
          ## Select a candidate for a new function of matching range type. This can of course result
          ## in a candidate function with a different domain type. If this happens the mutation is
          ## simply aborted, because searching again for a matching function would cost too much time...
          oldfunc <- expr[[1]]
          newfunc <- if (hasStype(oldfunc)) {
            oldfuncType <- sType(oldfunc)
            oldfuncRangeType <- rangeTypeOfType(oldfuncType)
            newfunccandidate <- toName(randelt(funcset$byRange[[oldfuncRangeType$string]],
                                               prob = attr(funcset$byRange[[oldfuncRangeType$string]], "probabilityWeight")))
            newfunccandidateType <- sType(newfunccandidate)
            if (identical(newfunccandidateType, oldfuncType)) newfunccandidate else oldfunc
          } else {
            newfunccandidate <- toName(randelt(funcset$all, prob = attr(funcset$all, "probabilityWeight")))
            if (arity(newfunccandidate) == arity(oldfunc)) newfunccandidate else oldfunc
          }
        } else expr[[1]]
        restExpr <- rest(expr)
        mutatedExpr <- as.call(append(mutatedLabel, Map(mutateExpressionChangeLabel, restExpr)))
        mutatedExpr
      }
    } else if (is.numeric(expr)) {
      if (currentNode %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
        mutatedExpr <- expr + rnorm(1)
        mutatedExpr
      } else expr
    } else if (is.logical(expr)) {
      if (currentNode %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
        mutatedExpr <- as.logical(rbinom(1, 1, 0.5))
        mutatedExpr
      } else expr
    } else stop("mutateChangeLabel: Unsupported expression: ", expr, ".")
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutateExpressionChangeLabel(body(func))
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateChangeLabel) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateInsertSubtree <- function(func, funcset, inset, conset,
                                strength = 1,
                                subtreeDepth = 2,
                                breedingFitness = function(individual) TRUE,
                                breedingTries = 50) {
  ## Subtrees are inserted by replacing a random leaf...
  numberOfLeaves <- funcLeaves(func)
  sampledMutationPoints <- sample.int(numberOfLeaves, replace = FALSE,
                                      size = min(strength, numberOfLeaves))
  currentLeaf <- 0 # here, local mutable state is more efficient than recursion
  mutateExpressionInsertSubtree <- function(expr) {
    if (is.call(expr)) {
      if (identical(expr[[1]], as.symbol("function"))) {
        stop("mutateInsertSubtree: Support for anonymous function nodes is not implemented.") # TODO
      } else {
        restExpr <- rest(expr)
        mutatedExpr <- as.call(append(expr[[1]], Map(mutateExpressionInsertSubtree, restExpr)))
        mutatedExpr
      }
    } else if (is.symbol(expr) || is.numeric(expr) || is.logical(expr)) {
      currentLeaf <<- currentLeaf + 1
      if (currentLeaf %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
        #message("INSERT SUBTREE at leaf#: ", currentLeaf, " -- ", expr) # DEBUG
        if (hasStype(expr)) {
          type <- sType(expr)
          newSubtree <- randexprTypedFull(type, funcset, inset, conset,
                                          maxdepth = subtreeDepth, constprob = 0.2)
          newSubtree
        } else {
          newSubtree <- randexprFull(funcset, inset, conset,
                                     maxdepth = subtreeDepth, constprob = 0.2)
          newSubtree
        }
      } else expr
    } else stop("mutateInsertSubtree: Unsupported expression: ", expr, ".")
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutateExpressionInsertSubtree(body(func))
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateInsertSubtree) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateDeleteSubtree <- function(func, funcset, inset, conset,
                                strength = 1,
                                subtreeDepth = 2,
                                constprob = 0.2,
                                breedingFitness = function(individual) TRUE,
                                breedingTries = 50) {
  ## Subtrees are deleted by deleting a random subtree of depth subtreeDepth...
  numberOfSubtreesOfMatchingDepth <- funcCount(func, function(node) exprDepth(node) == subtreeDepth)
  sampledMutationPoints <- sample.int(numberOfSubtreesOfMatchingDepth, replace = FALSE,
                                      size = min(strength, numberOfSubtreesOfMatchingDepth))
  currentMatchingSubtree <- 0 # here, local mutable state is more efficient than recursion
  mutateExpressionDeleteSubtree <- function(expr) {
    if (exprDepth(expr) == subtreeDepth) {
      currentMatchingSubtree <<- currentMatchingSubtree + 1
      if (currentMatchingSubtree %in% sampledMutationPoints && runif(1) > buildingBlockTag(expr)) {
        #message("DELETE SUBTREE at node#: ", currentMatchingSubtree, " -- ", expr) # DEBUG
        if (hasStype(expr)) {
          typeString <- rangeTypeOfType(sType(expr))$string
          newLeaf <- randterminalTyped(typeString, inset, conset, constprob)
          newLeaf
        } else {
          newLeaf <- if (runif(1) <= constprob) { # create constant
            constfactory <- randelt(conset$all, prob = attr(conset$all, "probabilityWeight"))
            constfactory()
          } else { # create input variable
            toName(randelt(inset$all, prob = attr(inset$all, "probabilityWeight")))
          }
          newLeaf
        }
      } else expr
    } else if (is.call(expr)) {
      if (identical(expr[[1]], as.symbol("function"))) {
        stop("mutateDeleteSubtree: Support for anonymous function nodes is not implemented.") # TODO
      } else {
        restExpr <- rest(expr)
        mutatedExpr <- as.call(append(expr[[1]], Map(mutateExpressionDeleteSubtree, restExpr)))
        mutatedExpr
      }
    } else if (is.symbol(expr) || is.numeric(expr) || is.logical(expr)) {
      expr
    } else stop("mutateDeleteSubtree: Unsupported expression: ", expr, ".")
  }
  doMutation <- function() {
    mutant <- new.function(envir = funcset$envir)
    formals(mutant) <- formals(func)
    body(mutant) <- mutateExpressionDeleteSubtree(body(func))
    mutant
  }
  breed(doMutation, breedingFitness, breedingTries)
}
class(mutateDeleteSubtree) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateChangeDeleteInsert <- function(func, funcset, inset, conset,
                                     strength = 1,
                                     subtreeDepth = 2,
                                     constprob = 0.2,
                                     iterations = 1,
                                     changeProbability = 1/3,
                                     deleteProbability = 1/3,
                                     insertProbability = 1/3,
                                     breedingFitness = function(individual) TRUE,
                                     breedingTries = 50) {
  stopifnot(iterations >= 1)
  if (iterations == 1) {
    mutateOperator <- randelt(c(mutateChangeLabel, mutateDeleteSubtree, mutateInsertSubtree),
                              prob = c(changeProbability, deleteProbability, insertProbability))
    do.call.ignore.unused.arguments(mutateOperator,
                                    list(func = func, funcset = funcset, inset = inset, conset = conset,
                                         strength = strength, subtreeDepth = subtreeDepth, constprob = constprob,
                                       breedingFitness = breedingFitness, breedingTries = breedingTries))
  } else {
    doMutate <- function(func) {
      mutateOperator <- randelt(c(mutateChangeLabel, mutateDeleteSubtree, mutateInsertSubtree),
                                prob = c(changeProbability, deleteProbability, insertProbability))
      do.call.ignore.unused.arguments(mutateOperator,
                                      list(func = func, funcset = funcset, inset = inset, conset = conset,
                                           strength = strength, subtreeDepth = subtreeDepth, constprob = constprob,
                                           breedingFitness = breedingFitness, breedingTries = breedingTries))
    }
    iterate(iterations, doMutate, func)
  }
}
class(mutateChangeDeleteInsert) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateDeleteInsert <- function(func, funcset, inset, conset,
                               strength = 1,
                               subtreeDepth = 2,
                               constprob = 0.2,
                               iterations = 1,
                               deleteProbability = 0.5,
                               insertProbability = 0.5,
                               breedingFitness = function(individual) TRUE,
                               breedingTries = 50)
  mutateChangeDeleteInsert(func, funcset, inset, conset,
                           strength = strength, subtreeDepth = subtreeDepth,
                           constprob = constprob, iterations = iterations,
                           changeProbability = 0,
                           deleteProbability = deleteProbability,
                           insertProbability = insertProbability,
                           breedingFitness = breedingFitness,
                           breedingTries = breedingTries)
class(mutateDeleteInsert) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateFuncFast <- function(funcbody, funcset, mutatefuncprob = 0.1)
  .Call("mutate_functions_R", funcbody, mutatefuncprob, funcset$all, as.integer(funcset$arities))
class(mutateFuncFast) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateSubtreeFast <- function(funcbody, funcset, inset, constmin, constmax, insertprob, deleteprob, subtreeprob, constprob, maxsubtreedepth)
  .Call("mutate_subtrees_R", funcbody, insertprob, deleteprob,
                             funcset$all, as.integer(funcset$arities),
                             inset$all, constmin, constmax, subtreeprob, constprob, as.integer(maxsubtreedepth))
class(mutateSubtreeFast) <- c("mutationOperator", "function")

##' @rdname expressionMutation
##' @export
mutateNumericConstFast <- function(funcbody, mutateconstprob = 0.1, mu = 0.0, sigma = 1.0)
  .Call("mutate_constants_normal_R", funcbody, mutateconstprob, mu, sigma)
class(mutateNumericConstFast) <- c("mutationOperator", "function")

