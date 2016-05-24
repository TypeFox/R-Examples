## GIES algorithm
##
## Author: Alain Hauser <alhauser@ethz.ch>
## $Id: gies.R 339 2015-07-22 11:25:06Z mmaechler $
###############################################################################

##################################################
## Auxiliary functions for simulations
##################################################

#' Randomly generates a Gaussian causal model >>>  ../man/r.gauss.pardag.Rd
#'
#' @param p number of vertices
#' @param prob probability of inserting an edge between two given
#'                    vertices
#' @param top.sort indicates whether the produced DAG should be
#'                    topologically sorted
#' @param normalize indicates whether weights and error variances
#'                    should be normalized s.t. the diagonal of the
#'                    corresponding covariance matrix is 1. Note that
#'                    weights and error variances can then lie outside
#'                    the boundaries specified below!
#' @param lbe lower bound of edge weights. Default: 0.1
#' @param ube upper bound of edge weights. Default: 1
#' @param neg.coef indicates whether also negative edge weights should
#'                    be sampled
#' @param labels
#' @param lbv lower bound of vertex variance. Default: 0.5
#' @param ubv upper bound of vertex variance. Default: 1
#' @return  an instance of gauss.pardag
r.gauss.pardag <- function(p,
    prob,
    top.sort = FALSE,
    normalize = FALSE,
    lbe = 0.1,
    ube = 1,
    neg.coef = TRUE,
    labels = as.character(1:p),
    lbv = 0.5,
    ubv = 1)
{
  ## Error checking
  stopifnot(is.numeric(p), length(p) == 1, p >= 2,
      is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
      is.numeric(lbe), is.numeric(ube), lbe <= ube,
      is.logical(neg.coef),
      is.numeric(lbv), is.numeric(ubv), lbv <= ubv,
      is.character(labels), length(labels) == p)

  ## Create list of nodes, edges and parameters
  edL <- as.list(labels)
  names(edL) <- labels

  ## Create list of parameters; first entry: error variances
  pars <- as.list(runif(p, min = lbv, max = ubv))
  names(pars) <- labels

  ## Create topological ordering
  top.ord <- if (top.sort) 1:p else sample.int(p)

  ## Sample edges and corresponding coefficients, respecting the generated
  ## topological ordering
  for (i in 2:p) {
    ii <- top.ord[i]
    parentCount <- rbinom(1, i - 1, prob)
    edL[[ii]] <- top.ord[sample.int(i - 1, size = parentCount)]
    weights <- runif(parentCount, min = lbe, max = ube)
    if (neg.coef)
      weights <- weights * sample(c(-1, 1), parentCount, replace = TRUE)
    pars[[ii]] <- c(pars[[ii]], 0, weights)
  }
  edL[[top.ord[1]]] <- integer(0)
  pars[[top.ord[1]]] <- c(pars[[top.ord[1]]], 0)

  ## Create new instance of gauss.pardag
  result <- new("GaussParDAG", nodes = labels, in.edges = edL, params = pars)

  ## Normalize if requested
  if (normalize) {
    H <- diag(result$cov.mat())
    result$set.err.var(result$err.var() / H)
    H <- sqrt(H)
    for (i in 1:p)
      if (length(edL[[i]]) > 0)
        result$.params[[i]][-c(1, 2)] <- pars[[i]][-c(1, 2)] * H[edL[[i]]] / H[i]
  }

  ## Validate and return object
  validObject(result)
  result
}

#' Simulates independent observational or interventional data for a
#' specified interventions from a Gaussian causal model
#'
#' @param   n         number of data samples
#' @param   object    an instance of gauss.pardag
#' @param   target    intervention target
#' @param   target.value    value of intervention targets
rmvnorm.ivent <- function(n, object, target = integer(0), target.value = numeric(0))
{
  p <- object$node.count()
  ## Error checking
  stopifnot(length(target) == 0 || (1 <= min(target) && max(target) <= p))
  stopifnot((is.vector(target.value) && length(target.value) == length(target)) ||
            (is.matrix(target.value) && dim(target.value) == c(n, length(target))))

  ## Simulate error terms
  sigma <- sqrt(object$err.var())
  mu <- object$intercept()
  Y <- matrix(rnorm(n*p, mu, sigma), nrow = p, ncol = n)

  ## Insert intervention values
  Y[target, ] <- target.value

  ## Calculate matrix of structural equation system
  A <- - t(object$weight.mat(target))
  diag(A) <- 1.

  ## Solve linear structural equations
  t(solve(A, Y))
}


##################################################
## Structure learning algorithms
##################################################
caus.inf <- function(algorithm, p, targets, score, ...)
{
  essgraph <- new("EssGraph", nodes = as.character(1:p), targets = targets, score = score)
  if (essgraph$caus.inf(algorithm, ...))
    list(essgraph = essgraph, repr = essgraph$repr())
}

gies <- function(p, targets, score, fixedGaps = NULL,
                 turning = TRUE, maxDegree = integer(0), verbose = FALSE, ...)
    caus.inf("GIES", p, targets, score, fixedGaps = fixedGaps,
             turning = turning, maxDegree = maxDegree, verbose = verbose, ...)

ges <- function(p, score, fixedGaps = NULL,
                turning = TRUE, maxDegree = integer(0), verbose = FALSE, ...)
    caus.inf("GIES", p, list(integer(0)), score, fixedGaps = fixedGaps,
             turning = turning, maxDegree = maxDegree, verbose = verbose, ...)

## TODO: make sure that the "representative" in the result is actually the last
## visited DAG instead of a random representative; adapt documentation accordingly
gds <- function(p, targets, score, verbose = FALSE, ...)
    caus.inf("GDS", p, targets, score, verbose = verbose, ...)

simy <- function(p, targets, score, verbose = FALSE, ...)
    caus.inf("SiMy", p, targets, score, verbose = verbose, ...)

#' Create a list of targets and a vector of target indices out of a
#' matrix indicating interventions
#'
#' @param 	A		a n x p boolean matrix; A[i, j] is TRUE iff vertex j is intervened
#' 							in data point i
#' @return 	list with two entries, "targets" and "target.index".
#' 					targets is a list of unique intervention targets
#' 					target.index is a vector of size n; the intervention target of data point
#' 					i is given by targets[[target.index[i]]].
mat2targets <- function(A)
{
  targets.raw <- as.list(apply(A, 1, which))
  targets <- unique(targets.raw)
  list(targets = targets, target.index = match(targets.raw, targets))
}

dag2essgraph <- function(dag, targets = list(integer(0))) {
  new("EssGraph",
      nodes = dag$.nodes,
      in.edges = .Call("dagToEssentialGraph", dag$.in.edges, targets),
      targets = targets)
}

#' Fast version of "gaussCItest", implemented in C++
gaussCItest.fast <- function(x, y, S, suffStat)
  .Call("condIndTestGauss", x, y, S, suffStat$n, suffStat$C)

