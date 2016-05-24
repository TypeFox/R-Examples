################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module B: rags2ridges Fused
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

##------------------------------------------------------------------------------
##
## (Hidden) Support functions
##
##------------------------------------------------------------------------------

isSymmetricPD <- function(M) {
  ##############################################################################
  # - Test if matrix is symmetric postive definite
  # - M > A numeric matrix
  #
  # NOTES:
  # - Returns TRUE if the matrix is symmetric positive definite and FALSE if not
  ##############################################################################

  nm <- deparse(substitute(M))
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(nm, " is not a numeric matrix")
  }
  if (!isSymmetric(M)) {
    stop(nm, " is not a symmetric matrix")
  }

  chM <- try(chol(M), silent = TRUE)  # M is P.D. iff it has a Cholesky decomp.
  if (inherits(chM, "try-error")) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}



isSymmetricPSD <- function(M, tol = 1e-4) {
  ##############################################################################
  # - Test if matrix is symmetric postive semi-definite
  # - M > A numeric matrix
  #
  # NOTES:
  # - Returns TRUE if the matrix is symmetric positive definite and FALSE if not
  # - In practice, it tests if all eigenvalues are larger than -tol*|l| where
  #   l is the largest eigenvalue.
  # - See http://scicomp.stackexchange.com/questions/12979/
  #          testing-if-a-matrix-is-positive-semi-definite
  ##############################################################################

  nm <- deparse(substitute(M))
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(nm, " is not a numeric matrix")
  }
  if (!isSymmetric(M)) {
    stop(nm, " is not a symmetric matrix")
  }

  evals <- eigen(M, symmetric = TRUE)$values # M is PSD iff eigenvals >= 0 SLOW!
  if (all(evals > -tol*abs(max(evals)))) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}



is.Xlist <- function(Xlist, Ylist = FALSE, semi = FALSE) {
  ##############################################################################
  # - Test if generic fused list arguments (such as Slist, Tlist, Plist)
  #   are properly formatted
  # - Xlist > A list of covariance matrices or matrices.
  # - Ylist > Is the supplied list a "Ylist"?
  # - semi  > Should the matrices be tested for postive (semi) definiteness
  #           If TRUE postive definiteness is tested.
  #
  # NOTES:
  # - Returns TRUE if all tests are passed, throws error if not.
  ##############################################################################

  xlist <- deparse(substitute(Xlist))
  if (!is.list(Xlist)) {
    stop(xlist, " should be a list")
  }
  if (!all(sapply(Xlist, is.matrix))) {
    stop("All elements of ", xlist, " should be matrices")
  }
  if (!all(sapply(Xlist, is.numeric))) {
    stop("All elements of ", xlist, " should be numeric matrices")
  }
  if (length(unique(c(sapply(Xlist, dim)))) != 1L) {
    stop("All matrices in ", xlist,
         " should be square and have the same size.")
  }
  if (semi) {
    if (!all(sapply(Xlist, isSymmetricPSD))) {
      stop("All matrices in ", xlist, " should be symmetric and positive ",
           "semi definite.")
    }
  } else {
    if (!all(sapply(Xlist, isSymmetricPD))) {
      stop("All matrices in ", xlist, " should be symmetric and positive ",
           "definite.")
    }
  }
  if (!all(sapply(seq_along(Xlist),
                  function(i) identical(dimnames(Xlist[[1]]),
                                        dimnames(Xlist[[i]]))))) {
    stop("dimnames of the elements of ", xlist, " are not identical")
  }

  return(TRUE)
}



default.target.fused <- function(Slist, ns, type = "DAIE", by, ...) {
  ##############################################################################
  # - Generate a list of (data-driven) targets to use in fused ridge estimation
  # - A nice wrapper for default.target
  # - Slist > A list of covariance matrices
  # - ns    > A numeric vector of sample sizes corresponding to Slist
  # - type  > A character giving the choice of target. See default.target.
  # - by    > A character or factor of the same length as Slist specifying
  #           which groups should share target. If omitted, each class
  #           has a unique target.
  # - ...   > Arguments passed to default.target.
  #           See also default.target
  ##############################################################################

  stopifnot(is.list(Slist))
  stopifnot(is.numeric(ns))
  stopifnot(length(Slist) == length(ns))

  if (missing(by)) {
    by <- seq_along(Slist)
  }
  stopifnot(length(by) == length(Slist))
  by <- as.character(by)
  stopifnot(!anyNA(by) && is.character(by))

  Tlist <- vector("list", length(Slist))
  names(Tlist) <- names(Slist)
  for (lvl in unique(by)) {
    get <- lvl == by
    pooled <- pooledS(Slist, ns, subset = get)
    Tpool <- default.target(pooled, type = type, ...)
    Tlist[get] <- replicate(sum(get), Tpool, simplify = FALSE)
  }

  return(Tlist)
}



createS <- function(n, p,
                    topology = "identity",
                    dataset = FALSE,
                    precision = FALSE,
                    nonzero = 0.25,
                    m = 1L,
                    banded.n = 2L,
                    invwishart = FALSE,
                    nu = p + 1,
                    Plist) {
  ##############################################################################
  # - Simulate one or more random symmetric square (or datasets) matrices from
  #   various models.
  # - n          > A vector of sample sizes
  # - p          > An integer giving the dimension. p should be greater than
  #                or equal to 2.
  # - topology   > character. Specify the topology to simulate data from.
  #                See details.
  # - dataset    > logical. Should dataset instead of its sample covariance be
  #                returned? Default is FALSE.
  # - precision  > logical. Should the constructed precision matrix be returned?
  # - nonzero    > numeric of length 1 giving the value of the non-zero entries
  #                for some topologies
  # - m          > integer. The number of conditionally independent subgraphs.
  #                I.e. the number of blocks.
  # - banded.n   > interger. The number of off-diagonal bands used if
  #                topology is "banded". Use as paremeter if topology is
  #                "Watt-Strogatz".
  # - invwishart > logical. If TRUE the constructed precision matrix is
  #                used as the scale matrix of an inverse Wishart distribution
  #                and class covariance matrices are drawn from this
  #                distribution.
  # - nu         > The degrees of freedom in the inverse wishart distribution.
  #                A large nu implies high class homogeneity.
  #                Must be greater than p + 1.
  # - Plist      > A list of user-supplied precision matrices. Should be the
  #                same length as n.
  #
  # NOTES:
  # - Returns a list of matrices if length(n) > 1. The output is simplified if
  #   n has length 1 where only the matrix is returned
  ##############################################################################

  if (missing(p) && !missing(Plist)) {
    p <- nrow(Plist[[1]])
  }
  stopifnot(p > 1)
  stopifnot(m >= 1)
  G <- length(n)

  if (dataset && precision) {
    stop("dataset and precision cannot be TRUE at the same time.")
  }

  if (invwishart && missing(nu)) {
    stop("argument 'nu' is missing. Supply the degrees of freedom 'nu' for ",
         "the inverse Wishart distribution.")
  }

  topology <- match.arg(topology,
                        c("identity", "star", "clique", "complete",
                          "chain", "banded", "Barabassi", "small-world",
                          "scale-free", "Watts-Strogatz", "random-graph",
                          "Erdos-Renyi"))

  # Construct the precision matrix "constructor"
  if (topology == "identity") {

    submat <- function(p) diag(p)

  } else if (topology == "star") {

    submat <- function(p) {
      subP <- diag(p)
      subP[1, seq(2, p)] <- subP[seq(2, p), 1] <- 1/seq(2, p)
      return(subP)
    }

  } else if (topology == "chain") {

    submat <- function(p) {
      s <- seq_len(p - 1)
      subP <- diag(p)
      subP[cbind(s, s + 1)] <- subP[cbind(s + 1, s)] <- nonzero
      return(subP)
    }

  } else if (topology == "clique" || topology == "complete") {

    submat <- function(p) {
      subP <- diag(p)
      subP[lower.tri(subP)] <- subP[upper.tri(subP)]  <- nonzero
      return(subP)
    }

  } else if (topology == "banded") {

    submat <- function(p) {
      if (banded.n > p) {
        stop("The number of bands cannot exceed the dimension of each block")
      }
      subP <- diag(p)
      for (j in seq(1, banded.n)) {
        s <- seq_len(p - j)
        subP[cbind(s, s + j)] <- subP[cbind(s + j, s)] <- 1/(j + 1)
      }
      return(subP)
    }

  } else if (topology == "Barabassi" || topology == "scale-free") {

    submat <- function(p) {
      G <- barabasi.game(p, power = 1, directed = FALSE)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else if (topology == "Watts-Strogatz" || topology == "small-world") {

    submat <- function(p) {
      G <- watts.strogatz.game(1, p, banded.n, 0.05)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else if (topology == "Erdos-Renyi" || topology == "random-graph") {

    submat <- function(p) {
      G <- erdos.renyi.game(p, 1/p)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else {

    stop("requested topology not implemented yet.")

  }

  # Construct the block split
  blocks <- split(seq_len(p), ceiling(m*seq_len(p)/p))

  # Fill in blocks to construct full precisions
  P <- diag(p)
  for (b in blocks) {
    P[b, b] <- submat(length(b))
  }
  if (rcond(P) < sqrt(.Machine$double.eps)) {  # Check condition number
    warning("The generated precision matrix has a very high condition number ",
            "and the generated data might be unreliable.")
  }
  S <- solve(P)

  # Construct names
  n.letters <- which.max(p <= 26^(1:3))
  x <- expand.grid(rep(list(LETTERS), n.letters))
  nms <- do.call(paste0, x)

  # Construct list to fill and iterate over all classes
  ans <- vector("list", G)
  names(ans) <- paste0("class", seq_len(G))
  for (g in seq_len(G)) {

    if (!missing(Plist)) {
      stopifnot(length(Plist) == length(n))
      stopifnot(nrow(Plist[[g]]) == ncol(Plist[[g]]))
      stopifnot(nrow(Plist[[g]]) == p)

      Sg <- solve(Plist[[g]])
    } else if (invwishart) {
      stopifnot(nu - p - 1 > 0)
      Sg <- drop(.armaRInvWishart(n = 1, psi = (nu - p - 1)*S, nu = nu))
    } else {
      Sg <- S
    }

    if (precision) {

      if (invwishart) {
        ans[[g]] <- solve(Sg)
      } else {
        ans[[g]] <- P
      }

    } else {

      ans[[g]] <- rmvnormal(n = n[g], mu = rep(0, p), sigma = Sg)

      if (!dataset) {
        ans[[g]] <- covML(ans[[g]])
      }

    }

    if (p <= 17576) {  # Only give names for "small" dimensions
      colnames(ans[[g]]) <- nms[1:p]
      if (!dataset) {
        rownames(ans[[g]]) <- nms[1:p]
      }
    }
  }

  if (G == 1) {  # Simplify output if ns is length 1
    ans <- ans[[1]]
  }

  return(ans)
}



getKEGGPathway <- function(kegg.id) {
  ##############################################################################
  # - Download information and graph object of a given kegg pathway.
  # - kegg.id  > The kegg id, e.g. "map04210", "map04064", "map04115". Can
  #              be prefixed with "path:".
  #
  # NOTES:
  # - The igraph objects can be obtained with igraph::igraph.from.graphNEL().
  # - The moral graph can be obtained with gRbase::moralize().
  # - To obtain the adjacency matrix, use gRbase::as.adjMAT() or
  #   igraph::get.adjacency()
  ##############################################################################

  if (!requireNamespace("KEGGgraph", quietly=TRUE)) {

    stop("This function requires the bioconductor package 'KEGGgraph' and its.",
         "\ndependencies. To use it, install KEGGgraph to use it by running:\n",
         'source("http://bioconductor.org/biocLite.R")\n',
         'biocLite("KEGGgraph")\n')
  }

  # Download
  tmp.file <- paste0(tempfile(), ".kgml")
  kgml <- KEGGgraph::retrieveKGML(gsub("path:", "", kegg.id), organism="hsa",
                                  destfile = tmp.file, method = "internal")

  # Pathway data.frame information
  df        <- KEGGgraph::parseKGML2DataFrame(tmp.file)
  df$to.e   <- KEGGgraph::translateKEGGID2GeneID(df$to)
  df$from.e <- KEGGgraph::translateKEGGID2GeneID(df$from)

  # Pathway igraph and graphNEL
  graph <- KEGGgraph::KEGGpathway2Graph(KEGGgraph::parseKGML(tmp.file))

  # Make sure that the graph is simple
  # A confirmed bug in KEGGgraph and now corrected in the development branch.
  # The following lines is a temporary work-around
  graph <- igraph.to.graphNEL(simplify(igraph.from.graphNEL(graph)))

  return(list(df = df, graph = graph))

}



.parents <- function(node, graph) {
  ##############################################################################
  # - Function for extracting the parents of a node
  # - node  > A characther giving the node name in the graph
  # - graph > The graph
  #
  # NOTES:
  # - Alternative to gRbase::parents()
  ##############################################################################

  if (!requireNamespace("graph")) {
    stop("package 'graph' needed for this function.")
  }

  is.child <- sapply(graph::edges(graph), function(n) node %in% n)
  return(graph::nodes(graph)[is.child])
}



kegg.target <- function(Y, kegg.id, method = "linreg", organism = "hsa",
                        graph = getKEGGPathway(kegg.id)$graph) {
  ##############################################################################
  # - Generate a target matrix from the KEGG database and pilot data.
  # - Requires a connection to the internet.
  # - Y        > The complete observation matrix of observations with variables
  #              in columns. The column names should be on the form e.g.
  #              "hsa:3988" ("<organism>:<Entrez id>"). It can however also be
  #              just the Entrez id with or without the post-fixed "_at" and
  #              then the specified organism will be assumed.
  # - kegg.id  > The kegg id, e.g. "map04210", "map04064", "map04115".
  # - method   > The method for estimating the non-zero entries moralized graph.
  #              Currently, only "linreg" is implemented.
  # - organism > A character giving the organism, default is
  #              "hsa" (homo-sapiens).
  # - graph    > A graphNEL object. Can be used to avoid repeatedly downloading
  #              the information.
  #
  # NOTES:
  # - See also default.target, and default.target.fused
  ##############################################################################

  method <- match.arg(method)
  stopifnot(length(organism) == 1L)

  if (!requireNamespace("KEGGgraph") && !requireNamespace("graph")) {

    stop("This function requires the bioconductor package 'KEGGgraph' and its.",
         "\ndependencies. To use it, install KEGGgraph to use it by running:\n",
         'source("http://bioconductor.org/biocLite.R")\n',
         'biocLite("KEGGgraph")\n')

  }

  # Check input

  correct.format <- grepl("^([[:alpha:]]+:)?[0-9]+(_at)?$", colnames(Y))
  s <- sum(!correct.format)
  if (s > 0) {
    wmsg <- paste("Found %d colnames of Y which incorrectly formatted.",
                  "They should be on the form <organism>:<entrez id> or",
                  "<entrez id> optionally be postfixed with '_at'")
    warning(sprintf(wmsg, s))
  }

  splt <- sapply(strsplit(colnames(Y), ":"),
                 function(x) if (length(x)==2) x[1] else organism)
  if (any(splt != organism)) {
    stop("The prefix does not always match the specified organism")
  }

  #
  # Download pathway and graphNEL object
  #

  stopifnot(is(graph, "graphNEL"))
  if (!is.dag(igraph.from.graphNEL(graph))) {
    warning("The graph obtained from KEGG is not acyclic. Results are only",
            " approximate.")
  }

  # Try to correct colnames (and save the old ones)
  colnames.org <- colnames(Y)
  colnames(Y) <- gsub("_at$", "", colnames(Y))  # Remove any _at post-fix,
  colnames(Y) <- gsub(paste0("^", organism, ":"), "", colnames(Y)) # rm prefix
  colnames(Y) <- paste0(organism, ":", colnames(Y)) # Put prefix back on all

  # Determine nodes/variables both on array and in pathway and subset
  common <- intersect(graph::nodes(graph), colnames(Y))
  ind <- match(common, colnames(Y))

  if (length(common) == 0) {
    stop("There were no gene IDs in pathway and supplied data. ",
         "Check that the column names are correctly formatted.")
  }
  g <- graph::subGraph(common, graph) #=removeNode(setdiff(nodes(G),common),G)
  Ysub <- Y[, common]

  # Center the data
  Ysub <- scale(Ysub, center = TRUE, scale = FALSE)

  if (method == "linreg") {

    # Intitialize precision matrix
    p <- graph::numNodes(g)
    prec <- matrix(0, p, p)
    rownames(prec) <- colnames(prec) <- common
    for (node in graph::nodes(g)) {

      pa.node <- .parents(node, g)
      # fit <- lm(Ysub[, node] ~ Ysub[, pa.node])  # Alternative computation
      # tausq <- summary(fit)$sigma^2
      fit <- lm.fit(x = cbind(Intercept = 1, Ysub[, pa.node, drop = FALSE]),
                    y = Ysub[, node])
      tausq <- (sum(fit$residuals^2)/(nrow(Ysub) - fit$rank))
      beta <- coef(fit)

      # Update entries [node, node]
      prec[node, node] <- prec[node, node] + 1/tausq

      # Update entries [node, pa(node)]
      prec[node, pa.node] <- prec[node, pa.node] - beta[pa.node]/tausq
      prec[pa.node, node] <- prec[node, pa.node]

      # Update entries [pa(node), pa(node)]
      prec[pa.node, pa.node] <-
        prec[pa.node, pa.node] + tcrossprod(beta[pa.node])/tausq

      if (anyNA(prec)) {
        stop("NAs were introduced")
      }
    }

    if (!isSymmetricPD(prec)) {
      warning("Constructed target is not symmetric PD")
    }

    rownames(prec) <- colnames(prec) <- colnames.org[ind]
    return(prec)

  } else {

    stop("No other methods are currently implmented")

  }

}



pooledS <- function(Slist, ns, subset = rep(TRUE, length(ns)), mle = TRUE) {
  ##############################################################################
  # - Computes the pooled covariance estimate
  # - Slist  > A list sample covariance matrices for each class
  # - ns     > A numeric vector of sample sizes of the same length as Slist.
  # - subset > logical vector the same length as Slist and ns giving the
  #            groups to pool over. Default is all.
  # - mle    > logical. If TRUE the biased MLE is used. If FALSE, the biased
  #            corrected estimate is used.
  ##############################################################################

  # Check input
  mle <- as.logical(mle)
  subset <- as.logical(subset)
  if (any(is.na(mle))) {
    stop("mle could not be coerced to a logical")
  }
  if (any(is.na(subset))) {
    stop("subset could not be coerced to a logical")
  }
  stopifnot(is.list(Slist) && length(Slist) == length(ns))
  stopifnot(is.logical(mle) && length(mle) == 1)
  stopifnot(is.logical(mle) && length(subset) == length(Slist))
  if (!any(subset)) {
    stop("argument subset must contain at least one TRUE entry.")
  }

  # Subsetting
  Slist <- Slist[subset]
  ns <- ns[subset]

  # Compute estimate
  ans <- .armaPooledS(Slist = Slist, ns = ns, mle = mle)
  dimnames(ans) <- dimnames(Slist[[1]])

  return(ans)
}



pooledP <- function(Plist, ns, subset = rep(TRUE, length(ns)), mle = TRUE) {
  ##############################################################################
  # - Computes the pooled precision estimate
  # - Plist  > A list (perhaps estimated) precision matrices for each class
  # - ns     > A numeric vector of sample sizes of the same length as Plist.
  # - subset > logical vector the same length as Plist and ns giving the
  #            groups to pool over. Default is all.
  # - mle    > logical. If TRUE the biased MLE is used. If FALSE, the biased
  #            corrected estimate is used.
  ##############################################################################

  # Check input
  mle <- as.logical(mle)
  subset <- as.logical(subset)
  if (any(is.na(mle))) {
    stop("mle could not be coerced to a logical")
  }
  if (any(is.na(subset))) {
    stop("subset could not be coerced to a logical")
  }
  stopifnot(is.list(Plist) && length(Plist) == length(ns))
  stopifnot(is.logical(mle) && length(mle) == 1)
  stopifnot(is.logical(mle) && length(subset) == length(Plist))
  if (!any(subset)) {
    stop("argument subset must contain at least one TRUE entry.")
  }

  # Subsetting
  Plist <- Plist[subset]
  ns <- ns[subset]

  # Compute estimate
  ans <- .armaPooledP(Plist = Plist, ns = ns, mle = mle)
  dimnames(ans) <- dimnames(Plist[[1]])

  return(ans)
}



KLdiv.fused <- function(MtestList, MrefList, StestList, SrefList, ns,
                        symmetric = FALSE) {
  ##############################################################################
  # - Function that calculates the "weighted" Kullback-Leibler divergence
  #   between multiple paired normal distributions
  # - MtestList > A list of mean vectors. Usually the empircal means.
  #               Assumed to be zero vectors if not supplied.
  # - MrefList  > A list of mean vectors. Usually the 'true'/reference/
  #               population means. Assumed to be zero vectors if not supplied.
  # - StestList > A list of (sample) covariance matrices.
  # - SrefList  > A list of (population) covariance matrices.
  # - ns        > A vector of sample sizes of the same length as Slist.
  #               Used as weights in the weighted mean.
  # - symmetric > logical indicating if original symmetric version of KL div.
  #               should be calculated.
  ##############################################################################

  if (missing(MtestList)) {
    MtestList <- replicate(length(StestList), rep(0, nrow(StestList[[1]])),
                           simplify = FALSE)
  }
  if (missing(MrefList)) {
    MrefList <- replicate(length(StestList), rep(0, nrow(StestList[[1]])),
                          simplify = FALSE)
  }
  KLdivs <- mapply(KLdiv, MtestList, MrefList, StestList, SrefList,
                   MoreArgs = list(symmetric = symmetric))

  return(sum(ns*KLdivs)/sum(ns))
}



.LL.fused <- function(Slist, Plist, ns){
  ##############################################################################
  # - Function that computes the value of the (negative) combined log-likelihood
  # - Slist > A list sample covariance matrices for each class
  # - Plist > A list of the same length as (Slist) of precision matrices
  #   (possibly regularized inverse of covariance or correlation matrices)
  # - ns > A vector of sample sizes of the same length as Slist.
  ##############################################################################

  LLs <- mapply(.LL, Slist, Plist)
  return(sum(ns*LLs))
}



.PLL.fused <- function(Slist, Plist, ns, Tlist, lambda){
  ##############################################################################
  # - Function that computes the value of the (negative) penalized combined
  #   log-likelihood
  # - Slist   > A list of G sample covariance matrices for each class
  # - Plist   > A list of G precision matrices with the same size as Slist.
  #             Possibly regularized inverse of covariance matrices.
  # - ns      > A vector of sample sizes of the same length as Slist.
  # - Tlist   > A list of G p.d. target matrices
  # - lambda  > A non-negative symmetric G by G penalty matrix
  ##############################################################################

  penalty <- 0
  for (g1 in seq_along(Slist)) {
    for (g2 in seq_len(g1)) {
      if (g1 == g2) { # Ridge penalty
        penalty <- penalty +
          lambda[g1, g1]*.FrobeniusLoss(Slist[[g1]], Tlist[[g1]])
      } else {  # Fused contribution
        penalty <- penalty +
          lambda[g1, g2]*.FrobeniusLoss(Slist[[g1]] - Tlist[[g1]],
                                        Slist[[g2]] - Tlist[[g2]])
      }
    }
  }
  penalty <- penalty/2

  ans <- .LL.fused(Slist, Plist, ns) + penalty
  return(ans)
}




##------------------------------------------------------------------------------
##
## Functions for the fused ridge estimator
##
##------------------------------------------------------------------------------

.init.ridgeP.fused <- function(Slist, ns, Tlist, lambda, ...) {
  ##############################################################################
  # - Internal function for selecting initial values for Plist
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A numeric non-negative symmetric G by G penalty matrix giving
  #             the penalties of the fused ridge estimator. The diagonal entries
  #             correspond to the class ridge penalites. The off-diagonal
  #             entries, lambda[g1, g2] say, determine the retainment of
  #             similarities between estimates in classes g1 and g2.
  #             If lambda is a single number, a diagonal penalty with lambda in
  #             the diagonal is used (lambda*diag(G)).
  #             If lambda is supplied as a numeric vector of two numbers,
  #             the first is used as a common ridge penalty and the second
  #             as a common fusion penalty.
  # - ...     > Arguments passed to .armaRidgeP
  ##############################################################################

#   Spool <- pooledS(Slist, ns)
#   if (length(unique(Tlist)) == 1L) { # If all targets equal
#
#     init.Plist <-
#       .armaRidgeP(Spool, target = Tlist[[1]], .trace(lambda)/sum(ns), ...)
#     init.Plist <- replicate(length(ns), init.Plist, simplify = FALSE)
#
#   } else {
#
#     init.Plist <- vector("list", length(ns))
#     for (i in seq_along(ns)) {
#       init.Plist[[i]] <-
#         .armaRidgeP(Spool, target = Tlist[[i]], lambda[i,i]/ns[i], ...)
#     }
#
#   }

  init.Plist <- default.target.fused(Slist, ns, type = "DAIE")

  names(init.Plist) <- names(Slist)
  return(init.Plist)
}



ridgeP.fused <- function(Slist,
                         ns,
                         Tlist = default.target.fused(Slist, ns),
                         lambda,
                         Plist,
                         maxit = 100L,
                         verbose = TRUE,
                         relative = TRUE,
                         eps = sqrt(.Machine$double.eps)) {
  ##############################################################################
  # - The user function for the fused ridge estimate for a given lambda.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A numeric non-negative symmetric G by G penalty matrix giving
  #             the penalties of the fused ridge estimator. The diagonal entries
  #             correspond to the class ridge penalites. The off-diagonal
  #             entries, lambda[g1, g2] say, determine the retainment of
  #             similarities between estimates in classes g1 and g2.
  #             If lambda is a single number, a diagonal penalty with lambda in
  #             the diagonal is used (lambda*diag(G)).
  #             If lambda is supplied as a numeric vector of two numbers,
  #             the first is used as a common ridge penalty and the second
  #             as a common fusion penalty.
  # - Plist   > A list of length G giving the initial estimates. If not supplied
  #             the ridge estimate of the pooled estimate is used.
  # - maxit   > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - relative> logical. Should the convergence criterion be on relative scale?
  # - eps     > numeric. A positive convergence criterion. Default is about 1e-8
  ##############################################################################

  stopifnot(length(Slist) == length(Tlist))
  G <- length(Slist)  # Number of groups

  # Hande special imputs of lambda
  if (!is.matrix(lambda) && is.numeric(lambda)) {
    if (length(lambda) == 1) {
      lambda <- diag(lambda, G)
    } else if (length(lambda) == 2) {
      tmp <- matrix(lambda[2], G, G)
      diag(tmp) <- lambda[1]
      lambda <- tmp
    } else {
      stop("If lambda is not a numeric matrix, it should have length 1 or 2.")
    }

  }
  if (!isSymmetric(lambda, check.attributes = FALSE)) {
    stop("lambda must be symmetric")
  }
  if (!all(lambda >= 0)) {
    stop("entries of lambda must be non-negative")
  }
  if (!all(diag(lambda) > 0)) {
    stop("the diagonal of lambda must be strictly postive")
  }

  # Initialize estimates with the regular ridges from the pooled covariance
  if (missing(Plist)) {
    Plist <- .init.ridgeP.fused(Slist, ns = ns, Tlist = Tlist, lambda = lambda)
  }
  stopifnot(length(Slist) == length(Plist))

  # Overwrite the starting estimate with the fused estimate
  Plist <- .armaRidgeP.fused(Slist = Slist, ns = ns, Tlist = Tlist,
                             lambda = lambda, Plist = Plist, maxit = maxit,
                             eps = eps, relative = relative, verbose = verbose)

  # Keep dimnames and names
  for (g in seq_along(Slist)) {
    dimnames(Plist[[g]]) <- dimnames(Slist[[g]])
  }
  names(Plist) <- names(Slist)

  return(Plist)
}




##------------------------------------------------------------------------------
##
## LOOCV (and approximation) for the fused setting
##
##------------------------------------------------------------------------------

.fcvl <- function(lambda, Ylist, Tlist, init.Plist, hotstart = FALSE, ...) {
  ##############################################################################
  # - (Internal) Computes the fused leave-one-out cross-validation loss for
  #   given penalty matrix
  # - lambda     > The G by G penalty matrix.
  # - Ylist      > A list of length G of matrices of observations with samples
  #                in the rows and variables in the columns. A least 2
  #                samples (rows) are needed in each entry.
  # - Tlist      > A list of length G of target matrices the same size
  #                as those of Plist. Default is given by default.target.
  # - init.Plist > Initial estimates used in .armaRidgeP.fused
  # - hotstart   > If TRUE, the latest estimate is used for each left out
  # - ...        > Arguments passed to .armaRidgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  Slist.org <- lapply(Ylist, covML)

  # If Plist is not supplied
  if (missing(init.Plist)) {
    init.Plist <- .init.ridgeP.fused(Slist.org, ns.org, Tlist, lambda)
  }

  slh <- numeric(sum(ns.org))  # To store LOOCV losses for each sample
  j <- 1
  for (g in seq_len(G)) {
    ns <- ns.org        # "Reset" number of samples in each group
    ns[g] <- ns[g] - 1  # Update sample size in g'th group
    this.Slist <- Slist.org
    for (i in seq_len(ns.org[g])) {
      this.Slist[[g]] <-
        covML(Ylist[[g]][-i, , drop = FALSE])
        # crossprod(Ylist[[g]][-i, , drop = FALSE])/ns[g]

      this.Plist <- .armaRidgeP.fused(Slist = this.Slist, ns = ns,
                                      Tlist = Tlist, lambda = lambda,
                                      Plist = init.Plist, verbose = FALSE, ...)

      Sig <- crossprod(Ylist[[g]][i,  , drop = FALSE])
      slh[j] <- ns[g]*.LL(Sig, this.Plist[[g]])

      if (hotstart) {
        init.Plist <- this.Plist
      }

      j <- j + 1
    }
  }

  return(mean(slh))
}



.kfcvl <- function(lambda, Ylist, Tlist, init.Plist, k, ...) {
  ##############################################################################
  # - (Internal) Computes the k-fold fused cross-validation loss for a penalty
  #   matrix. The data for each class is divided into k parts. The first part
  #   in each class is left out, the fused estimate is computed based on the
  #   remaning, and the loss is computed. Then this is repeated for the remaning
  #   parts.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > Initial estimates
  # - k       > The fold of the cross validation. I.e. k is the number of
  #             roughly equally sized parts the samples for each class are
  #             partitioned into.
  # - ...     > Arguments passed to .armaRidgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  if (min(ns.org) < k) {
    stop("The least class sample size is less than the specified k = ", k)
  }

  # If Plist is not supplied
  if (missing(init.Plist)) {
    Slist.org  <- lapply(Ylist, covML)
    init.Plist <- .init.ridgeP.fused(Slist.org, ns.org, Tlist, lambda)
  }

  # Split each class into k equally sized parts
  parts <- lapply(ns.org, function(n) sample(ceiling(k*seq_len(n)/n)))

  # Run through all k splits in each class
  slh <- matrix(0, k, G)
  for (i in seq_len(k)) {
    # Pick out ALL BUT the i'th fold in each class and compute estimate
    Ylist.i <- mapply(function(x, ind) x[ind != i, , drop = FALSE],
                      Ylist, parts, SIMPLIFY = FALSE)
    ns.i    <- sapply(Ylist.i, nrow)
    Slist.i <- lapply(Ylist.i, covML)
    Plist.i <- .armaRidgeP.fused(Slist = Slist.i, ns = ns.i, Tlist = Tlist,
                                 lambda = lambda, Plist = init.Plist,
                                 verbose = FALSE, ...)

    # Evaluate estimate with left out data:
    for (g in seq_len(G)) {
      Ylist.ig <- Ylist[[g]][parts[[g]] == i, , drop = FALSE]
      nig <- nrow(Ylist.ig)
      Sig <- crossprod(Ylist.ig)/nig
      slh[i, g] <- .LL(Sig, Plist.i[[g]])
    }

  }

  return(mean(slh))
}



.sfcvl <- function(lambda, Ylist, Tlist, Plist, ...) {
  ##############################################################################
  # - (Internal) Computes the "special" LOOCV loss for given penalty parameters
  # - Only updates the class estimate in which the sample is left out.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > Initial estimates
  # - ...     > Arguments passed to .armaRidgeP.fused and ridgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  Slist.org <- lapply(Ylist, covML)

  # If Plist is not supplied
  if (missing(Plist)) {
    Plist.org <-  ridgeP.fused(Slist = Slist.org, ns = ns.org, Tlist = Tlist,
                               lambda = lambda, verbose = FALSE, ...)
  } else {
    Plist.org <- Plist
  }

  slh <- numeric(sum(ns.org))
  j <- 1
  for (g in seq_len(G)) {
    ns <- ns.org        # "Reset" number of samples in each group
    ns[g] <- ns[g] - 1  # Update sample size in g'th group

    for (i in seq_len(ns.org[g])) {
      Plist <- Plist.org
      Slist <- Slist.org
      Slist[[g]] <- covML(Ylist[[g]][-i, , drop = FALSE])

      # Update only the estimate in group "g".
      # Note these exported C++ functions are index from g = 0
      if (sum(lambda) < 1e50) {
        Plist[[g]] <- .armaFusedUpdateI(g0 = g - 1,  Plist = Plist,
                                        Slist = Slist, Tlist = Tlist, ns = ns,
                                        lambda = lambda)
      } else {
        Plist[[g]] <- .armaFusedUpdateIII(g0 = g - 1,  Plist = Plist,
                                          Slist = Slist, Tlist = Tlist, ns = ns,
                                          lambda = lambda)
      }

      Sig <- crossprod(Ylist[[g]][i,  , drop = FALSE])
      slh[j] <- .LL(Sig, Plist[[g]])
      j <- j +1
    }
  }

  return(mean(slh))
}



.afcvl <- function(lambda, Ylist, Tlist, Plist, ...) {
  ##############################################################################
  # - (Internal) Computes the approximate LOOCV loss for at given penalty
  #   parameters.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > A list of inital values for the parameter estimates
  # - ...     > Arguments passed to ridgeP.fused
  ##############################################################################

  ns <- sapply(Ylist, nrow)
  G <- length(ns)
  Slist <- lapply(Ylist, covML)
  Plist <- ridgeP.fused(Slist = Slist, Tlist = Tlist, ns = ns,
                       lambda = lambda, verbose = FALSE, ...)
  n.tot <- sum(ns)
  nll <- .LL.fused(Slist = Slist, Plist = Plist, ns)/n.tot
  p <- nrow(Slist[[1]])
  bias <- 0

  # Implementation 1
  for (g in seq_along(ns)) {
    for (i in seq_len(ns[g])) {
      Sig <- crossprod(Ylist[[g]][i, , drop = FALSE])
      fac1 <- diag(nrow(Sig)) - Sig %*% Plist[[g]]
      fac2 <- Plist[[g]] %*% (Slist[[g]] - Sig) %*% Plist[[g]]
      bias <- bias  + sum(fac1 * fac2)/(2*n.tot)
    }
  }

#   # Implementation 2  (SVANTE)
#   for (g in seq_along(ns)) {
#     lambdabar <- (lambda + sum(lambdaF[g,-g]))/ns[g]
#     P1 <- Plist[[g]]
#     P2 <- P1 %*% P1
#     P2mP1 <- P2 - Plist[[g]]
#     P4mP3 <- P2mP1 %*% P2
#     for (i in seq_len(ns[g])) {
#       yig <- Ylist[[g]][i, , drop = TRUE]
#       b1 <- (yig %*% P2mP1) %*% yig
#       b2 <- lambdabar*((yig %*% P4mP3) %*% yig)
#       bias <- bias + (b1 + b2)
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#
#   # Implementation 3 (VUCACIC)
#   for (g in seq_along(ns)) {
#     for (i in seq_len(ns[g])) {
#       Sig <- crossprod(Ylist[[g]][i, , drop = FALSE])
#       bias <- bias +
#         sum((solve(Plist[[g]]) - Sig)*Plist[[g]]*(Slist[[g]]-Sig)*Plist[[g]])
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#   # Implementation 4
#   qf <- function(x, A) return(colSums(x * (A %*% x)))
#   for (g in seq_along(ns)) {
#     ng <- ns[g]
#     t1 <- p*log((ng-1)/ng)
#     for (i in seq_len(ng)) {
#       yik <- Ylist[[g]][i, ]
#       yOy <- qf(yik, Plist[[g]])
#       oneMinusyOy <- 1 - yOy/ng
#       t1 - log(oneMinusyOy) + (yOy^2/oneMinusyOy - yOy)/ng
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#
#   # Implementation 5
#   bias <- 0
#   for (g in seq_along(ns)) {
#     bias <- bias + ns[g]*sum(Plist[[g]]*(solve(Plist[[g]]) - Slist[[g]]))
#     bias <- bias + (1/rcond(Plist[[g]]))
#   }
#   bias <- bias/(2*n.tot)

  return(nll + bias)
}



.parseLambda <- function(lambda) {
  ##############################################################################
  # - A function to parse a character matrix that defines the class of penalty
  #   graphs and unique parameters for cross validation.
  # - Returns a data.frame of different indices for each level to be penalized
  #   equally.
  # - This data.frame is to be used to construct numeric matrices of penalties.
  # - lambda > A symmetric G by G character matrix defining the class of penalty
  #            matrices to cross validate over.
  #            Entries with NA, "" (the empty string), or "0" are
  #            interpreted as that that pair should omitted.
  #            Entries coercible to numeric are (in turn) interpreted as fixed
  ##############################################################################

  stopifnot(is.character(lambda))
  stopifnot(is.matrix(lambda))
  stopifnot(nrow(lambda) == ncol(lambda))

  # Handle special values
  lambda[is.na(lambda)] <- "0"
  lambda[lambda %in% ""] <- "0"

  parsedLambda <-
    data.frame(name = as.character(lambda), row = as.integer(row(lambda)),
               col = as.integer(col(lambda)), stringsAsFactors = FALSE)
  parsedLambda$val <- suppressWarnings({as.numeric(parsedLambda$name)})
  parsedLambda$fixed <- !is.na(parsedLambda$val)
  parsedLambda$variable <- !parsedLambda$fixed
  nf <- !parsedLambda$fixed
  parsedLambda$index[nf] <- as.numeric(as.factor(parsedLambda$name[nf]))

  u <- unique(subset(parsedLambda, select = -c(row, col)))

  attr(parsedLambda, "n.classes") <- nrow(lambda)
  attr(parsedLambda, "n.fixed") <- sum(u$fixed)
  attr(parsedLambda, "n.variables") <- nrow(u) - attr(parsedLambda, "n.fixed")
  return(parsedLambda)
}



.reconstructLambda <- function(lambdas, parsedLambda) {
  ##############################################################################
  # - Reconstruct the numeric penalty matrix lambda from a vector (lambdas)
  #   of penalties using the .parseLambda output.
  # - lambdas      > A numeric vector of the penalties. The length of lambdas
  #                  is the number of non-fixed entries in parsedLambda
  # - parsedLambda > A data.frame describing the penalty matrix.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  if (length(lambdas) != attributes(parsedLambda)$n.variables) {
    stop("The number of lambdas does not correspond with the number of",
         " non-fixed penalties given in parsedLambda")
  }

  G <- attributes(parsedLambda)$n.classes
  var <- parsedLambda$variable
  vals <- parsedLambda$val
  vals[var] <- lambdas[parsedLambda$index[var]]
  lambda <- matrix(vals, G, G)
  return(lambda)
}



.lambdasFromMatrix <- function(lambda.init, parsedLambda) {
  ##############################################################################
  # - Create the "lambdas" vector used in the optimizers from a numeric matrix.
  # - lambda.init  > A numeric matrix of the initial penalty matrix.
  # - parsedLambda > A data.frame describing the penalty matrix.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  u <- parsedLambda[!duplicated(parsedLambda$index), ]
  u <- u[!u$fixed, ]
  lambdas <- numeric(attributes(parsedLambda)$n.variables)
  stopifnot(length(lambdas) == nrow(u))
  lambdas[u$index] <- lambda.init[as.matrix(subset(u, select = c(row, col)))]

  # Try to reconstruct the given matrix
  relambda.init <- .reconstructLambda(lambdas, parsedLambda)
  if (!all(relambda.init == lambda.init)) {
    warning("The fixed penalties do not agree with the specified initial ",
            "penalty matrix.")
  }
  return(lambdas)
}



optPenalty.fused.grid <-
  function(Ylist, Tlist,
           lambdas = 10^seq(-5, 5, length.out = 15),
           lambdaFs = lambdas,
           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
           k = 10,
           verbose = TRUE,
           ...) {
  ##############################################################################
  # - Cross validation for the fused ridge estimator on a grid to determine
  #   optimal lambda and lambdaF.
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambdas     > A vector of ridge penalties
  # - lambdaFs    > A vector of fusion penalties
  # - cv.method   > The LOOCV type to use. Allowed values are LOOCV, aLOOCV,
  #                 sLOOCV, kCV for leave-one-out cross validation (LOOCV),
  #                 appproximate LOOCV, special LOOCV, and k-fold CV, resp.
  # - k           > Number of parts in k-fold CV. Only use if method is "kCV".
  # - ...         > Arguments passed to ridgeP.fused
  # - verbose     > logical. Print extra information. Defaults is TRUE.
  #
  # NOTES:
  # - The complete penalty graph is assumed (i.e. all ridge penalties
  #   equal and all fusion penalties equal)
  ##############################################################################

  cv.method <- match.arg(cv.method)

  if (missing(Tlist)) {  # If Tlist is not provided
    Tlist <- lapply(Ylist, function(Y) default.target(covML(Y)))
  }

  G <- length(Ylist)

  stopifnot(all(lambdas > 0))
  stopifnot(all(lambdaFs >= 0))
  stopifnot(all(is.finite(lambdas)))
  stopifnot(all(is.finite(lambdaFs)))

  if (cv.method == "LOOCV") {
    cvfunc <- .fcvl
  } else if (cv.method == "aLOOCV") {
    cvfunc <- .afcvl
  } else if (cv.method == "sLOOCV") {
    cvfunc <- .sfcvl
  } else if (cv.method == "kCV") {
    cvfunc <-  function(lambda, Ylist = Ylist, Tlist = Tlist, ...) {
      .kfcvl(lambda, Ylist = Ylist, Tlist = Tlist, k = k, ...)
    }
  } else {
    stop("cv.method not implmented.")
  }

  # Calculate CV scores
  if (verbose) {
    cat("Calculating cross-validated negative log-likelihoods...\n")
  }

  slh <- matrix(NA, length(lambdas), length(lambdaFs))
  dimnames(slh) <- list("lambdas" = lambdas, "lambdaFs" = lambdaFs)
  for (l1 in seq_along(lambdas)) {
    for (l2 in seq_along(lambdaFs)) {
      # Create penalty matrix
      lambda <- matrix(lambdaFs[l2], G, G)
      diag(lambda) <- lambdas[l1]

      # Evaluate loss
      slh[l1, l2] <- cvfunc(lambda = lambda, Ylist = Ylist, Tlist = Tlist, ...)

      if (verbose) {
        cat(sprintf("lambda = %.3e (%d), lambdaF = %.3e (%d), fcvl = %.3e\n",
                   lambdas[l1],  l1, lambdaFs[l2], l2, slh[l1, l2]))
      }
    }
  }

  output <- list(lambda = lambdas, lambdaF = lambdaFs, fcvl = slh)
  class(output) <- "optPenaltyFusedGrid"
  return(output)
}



print.optPenaltyFusedGrid <- function(x, ...) {
  with(x, print(fcvl))
  return(invisible(x))
}



plot.optPenaltyFusedGrid <- function(x, add.text = TRUE, add.contour = TRUE,
                                     col = rainbow(100, end = 0.8), ...) {
  with(x, {
    image(log(lambda), log(lambdaF), fcvl, col = col)
    if (add.contour) {
      contour(log(lambda), log(lambdaF), log(fcvl - min(fcvl) + 1), add = TRUE,
              nlevels = 15, col = "White", drawlabels = FALSE)
    }
    cols <- with(x, ifelse(fcvl == min(fcvl), "red", "black"))
    if (add.text) {
      text(log(lambda)[c(row(fcvl))], log(lambdaF)[c(col(fcvl))],
           sprintf("%0.1f", fcvl), cex = 0.7, col = cols)
    }
  })
  return(invisible(x))
}



optPenalty.fused.auto <-
  function(Ylist,
           Tlist,
           lambda,
           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
           k = 10,
           verbose = TRUE,
           lambda.init,
           maxit.ridgeP.fused = 1000,
           optimizer = "optim",
           maxit.optimizer = 1000,
           debug = FALSE,
           optim.control = list(trace = verbose, maxit = maxit.optimizer),
           ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   leave-one-out cross-validation using multi-dimensional optimization
  #   routines.
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambda      > A symmetric G by G character matrix defining the class of
  #                 penalty matrices to cross validate over. The unique elements
  #                 of lambda specify the penalties to determine. Pairs can be
  #                 left out using either of "", NA or "0". Penalties
  #                 can be fixed if they are coercible to numeric values,
  #                 e.g. "2.5".
  # - cv.method   > The LOOCV type to use. Allowed values are LOOCV, aLOOCV,
  #                 sLOOCV, kCV for leave-one-out cross validation (LOOCV),
  #                 appproximate LOOCV, special LOOCV, and k-fold CV, resp.
  # - k           > Number of parts in k-fold CV. Only use if method is "kCV".
  # - verbose     > logical. Should extra info be printed? Defaults to TRUE.
  # - lambda.init > A numeric penalty matrix of initial values.
  #                 If omitted, the function selects a starting values using
  #                 a common ridge penaltiy (determiend by 1D otimization)
  #                 setting all fusion penalties to zero.
  # - maxit.ridgeP.fused > integer. Max. number of iterations for ridgeP.fused
  # - optimizer          > character giving the stadard optimizer.
  #                        Either "optim" or "nlm".
  # - maxit.optimizer    > integer. Max. number of iterations for the optimizer.
  # - debug              > logical. If TRUE the raw output from the optimizer is
  #                        appended as an attribute to the output.
  # - optim.control      > Control arguments for optim.
  # - ...                > arguments passed to the optimizer.
  ##############################################################################

  cv.method <- match.arg(cv.method)
  G <- length(Ylist)

  if (missing(lambda)) {  # Handle missing lambda
    lambda <- matrix("fusion", G, G)
    diag(lambda) <- "ridge"
  }

  # Interpret given lambda
  parsedLambda <- .parseLambda(lambda)

  if (verbose) {  # Report interpretation
    n.fixed     <- attributes(parsedLambda)$n.fixed
    n.variables <- attributes(parsedLambda)$n.variables
    nonfix <- with(parsedLambda, paste(unique(name[!fixed]), collapse = ", "))
    fixed  <- with(parsedLambda, paste(unique(name[ fixed]), collapse = ", "))
    message("Found ", n.fixed + n.variables, " unique penalties of which ",
            n.fixed, " are interpreted as fixed and ", n.variables,
            " are to be determined by ", cv.method, ".\n",
            "Non-fixed parameters: ", nonfix,
            "\nFixed parameters: ", ifelse(n.fixed, fixed, "<none>"))
  }

  # Determine what loss function to use
  if (cv.method == "LOOCV") {
    cvfunc <- .fcvl
  } else if (cv.method == "aLOOCV") {
    cvfunc <- .afcvl
  } else if (cv.method == "sLOOCV") {
    cvfunc <- .sfcvl
  } else if (cv.method == "kCV") {
    cvfunc <-  function(lambda, Ylist = Ylist, Tlist = Tlist, ...) {
      .kfcvl(lambda, Ylist = Ylist, Tlist = Tlist, k = k, ...)
    }
  } else {
    stop("cv.method not implmented.")
  }

  # Construct optim/nlm objective function, parameters are assumed on log-scale
  cvl <- function(lambdas, ...) {
    elambdas <- exp(lambdas)
    lambda <- .reconstructLambda(elambdas, parsedLambda)
    return(cvfunc(lambda = lambda, Ylist = Ylist, Tlist = Tlist,
                  maxit = maxit.ridgeP.fused, ...))
  }

  if (missing(lambda.init)) {
    # Get somewhat sensible starting value for non-fixed diagonal entries
    # (ridge penalties) and by choosing off-diag lambda to be zero.
    f <- function(x) {
      lambdas <- suppressWarnings({.lambdasFromMatrix(diag(exp(x), G),
                                                      parsedLambda)})
      return(cvl(lambdas))
    }
    st <- optimize(f, lower = -20, upper = 20)$minimum
    lambdas.st <- suppressWarnings({log(.lambdasFromMatrix(diag(exp(st), G),
                                                           parsedLambda))})
    lambdas.st[lambdas.st == -Inf] <- -40
    lambdas.st[lambdas.st ==  Inf] <-  40

  } else {

    if (is.matrix(lambda.init) && is.numeric(lambda.init) &&
        isSymmetric(lambda.init)) {
      lambdas.st <- log(.lambdasFromMatrix(lambda.init, parsedLambda))
      lambdas.st[lambdas.st == -Inf] <- -40
      lambdas.st[lambdas.st ==  Inf] <-  40

    } else {
      stop("The supplied inital parameters must be a symmetric numeric matrix")
    }
  }

  if (optimizer == "optim") {

    ans <- optim(lambdas.st, fn = cvl, ..., control = optim.control)
    par <- ans$par
    val <- ans$value

    args <- list(...)
    if (!is.null(args$method)) {
      if (args$method != "SANN" && ans$convergence == 1) {
        warning("Iteration limit of optim had been reached.")
      }
      if (args$method == "Nelder-Mead" && ans$convergence == 10) {
        warning("Degeneracy of the Nelder-Mead simplex.")
      }
    }

  } else if (optimizer == "nlm") {

    ans <- nlm(cvl, lambdas.st, iterlim = maxit.optimizer, ...)
    par <- ans$estimate
    val <- ans$minimum

  }

  # Format optimal values
  opt.lambdas <- exp(par)
  opt.lambda <- .reconstructLambda(opt.lambdas, parsedLambda)
  dimnames(opt.lambda) <- dimnames(lambda)

  # Construct output
  res <- list(Plist = NA,
              lambda = opt.lambda,
              lambda.unique = NA,
              value = val)

  # Compute estimate at optimal values
  res$Plist <- ridgeP.fused(Slist = lapply(Ylist, covML),
                            ns = sapply(Ylist, nrow),
                            Tlist = Tlist, lambda = res$lambda,
                            maxit = maxit.ridgeP.fused, verbose = FALSE)

  res$lambda.unique <- unique(opt.lambda[lower.tri(opt.lambda, diag = TRUE)])
  names(res$lambda.unique) <- lambda[match(res$lambda.unique, opt.lambda)]

  if (debug) {
    attr(res, "optim.debug") <- ans
  }

  return(res)
}



optPenalty.fused <- function(Ylist, Tlist, lambda = default.penalty(Ylist),
                             cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
                             k = 10, grid = FALSE, ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   LOOCV using multi-dimensional optimization routines. A simple wrapper for
  #   optPenalty.fused.auto and optPenalty.fused.grid
  # - Ylist  > A list of length G of matrices of observations with samples
  #            in the rows and variables in the columns.
  # - Tlist  > A list of length G of target matrices the same size
  #            as those of Plist. Default is given by default.target.
  # - lambda > A G by G symmetric character matrix defining the class of penalty
  #            matrices to use. The unique elements of lambda specify the
  #            penalties to determine. Penalties can be fixed by using
  #            entries coercible to a numeric, e.g. "2.1".
  #            Pairs can be left out using either
  #            of "", NA or "0".
  # - method > The LOOCV type to use. Allowed values are LOOCV, aLOOCV,
  #            sLOOCV, kCV for leave-one-out cross validation (LOOCV),
  #            appproximate LOOCV, special LOOCV, and k-fold CV, resp.
  # - k      > Number of parts in k-fold CV. Only use if method is "kCV".
  # - grid   > logical Should grid based search be used? Default is FALSE.
  # - ...    > arguments passed to optPenalty.fused.auto and
  #            optPenalty.fused.grid
  ##############################################################################

  cv.method <- match.arg(cv.method)

  if (grid) {
    res <- optPenalty.fused.grid(Ylist = Ylist, Tlist = Tlist,
                                 cv.method = cv.method, k = k, ... )
  } else {
    res <- optPenalty.fused.auto(Ylist = Ylist, Tlist = Tlist,
                                 lambda = lambda,
                                 cv.method = cv.method, k = k,...)
  }
  return(res)
}




##------------------------------------------------------------------------------
##
## Automatic penalty matrix constructor
##
##------------------------------------------------------------------------------

.charAdjMat <- function(fac, name = "X", ordered = is.ordered(fac)) {
  ##############################################################################
  # - Create a complete character adjacency matrix from a factor. This function
  #   is used in the constructing the character penalty matrix in
  #   default.penalty.
  # - fac     > A factor of some length. Can be ordered.
  # - name    > A character giving the text which should appear in the adjacent
  #             entries. If not a character, the object name of fac is used.
  # - ordered > logical specifiying if fac should be interpreted as ordered.
  #
  # Examples:
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "Y")
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = NULL)
  #  rags2ridges:::.charAdjMat(ordered(factor(LETTERS[1:5])))
  ##############################################################################

  G <- nlevels(fac)
  if (is.character(name)) {
    lab <- name
  } else {
    lab <- deparse(substitute(fac))
  }
  if (ordered) {
    M <- matrix("", G, G)
    M[row(M) == (col(M) + 1)] <- lab
    M[row(M) == (col(M) - 1)] <- lab
  } else {
    M <- matrix(lab, G, G)
    diag(M) <- ""
  }
  rownames(M) <- colnames(M) <- levels(fac)
  return(M)
}



.char2num <- function(X) {
  ##############################################################################
  # - Create a character adjacency matrix to a numeric one
  # - X  > A character matrix where "" signify non-adjacency.
  #
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  # rags2ridges:::.char2num(A)
  ##############################################################################

  X[X != ""] <- "1"
  X[X == ""] <- "0"
  Y <- structure(as.numeric(X), dim = dim(X), dimnames = dimnames(X))
  return(Y)
}



.cartesianProd <- function(A, B) {
  ##############################################################################
  # - Construct the Cartesian product graph from two "character" matrices.
  # - A > A character matrix where "" signify non-adjacency.
  # - B > A character matrix where "" signify non-adjacency.
  #
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "X")
  # B <- rags2ridges:::.charAdjMat(factor(letters[4:5]), name = "Y")
  # rags2ridges:::.cartesianProd(A, B)
  ##############################################################################

  AI <- kronecker(.char2num(A), diag(nrow(B)), make.dimnames = TRUE)
  IB <- kronecker(diag(nrow(A)), .char2num(B), make.dimnames = TRUE)
  gprod <- AI + IB  # Kronecker sum

  ans <- kronecker(A, B, FUN = paste0, make.dimnames = TRUE)
  ans[!as.logical(gprod)] <- ""
  return(ans)
}



.tensorProd <- function(A, B) {
  ##############################################################################
  # - Construct the Tensor (or categorical) product graph from two "character"
  #   matrices.
  # - A > A character matrix where "" signify non-adjacency.
  # - B > A character matrix where "" signify non-adjacency.
  #
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "X")
  # B <- rags2ridges:::.charAdjMat(factor(letters[4:5]), name = "Y")
  # rags2ridges:::.tensorProd(A, B)
  ##############################################################################

  gprod <- kronecker(.char2num(A), .char2num(B), make.dimnames = TRUE)

  ans <- kronecker(A, B, FUN = paste0, make.dimnames = TRUE)
  ans[!as.logical(gprod)] <- ""
  return(ans)
}



default.penalty <- function(G, df,
                            type = c("Complete", "CartesianEqual",
                                     "CartesianUnequal", "TensorProd")) {
  ##############################################################################
  # - Select a one of standard penalty matrix types from a dataframe
  # - G     > The number of classes. Can also be list of length G such as
  #           the usual argument "Slist".
  #           Can be omitted if 'df' is given.
  # - df    > A data.frame with G rows and with factors in the columns.
  #           Columns of type character are coerced to factors.
  #           Can be omitted when 'type == "Complete"'.
  # - type  > A character giving the type of fused penalty graph to construct.
  #           Should be one of 'Complete' (default), 'CartesianEqual',
  #           'CartesianUnequal', or 'TensorProd' or an unique abbreviation
  #           hereof.
  #
  # NOTES:
  # - default.penalty can use ordered factors of df
  # - Setting type == 'Complete' is the complete penalty graph with equal
  #   penalties.
  # - Setting type == 'CartesianEqual' corresponds to a penalizing along each
  #   "direction" of factors with a common penalty.
  # - Setting type == 'CartesianUnequal' corresponds to a penalizing each
  #   direction of factors with individual penalties.
  ##############################################################################

  type <- match.arg(type)

  if (missing(G) && !missing(df)) {
    G <- nrow(df)
  }

  if (is.data.frame(G)) {
    df <- G
    G <- nrow(df)
  }

  if (missing(G) && !missing(df)) {
    G <- nrow(df)
  }

  if (is.list(G)) {
    G <- length(G)
  }

  if (missing(df)) {
    if (type != "Complete") {
      warning("No data.frame 'df' given and 'type' does not equal 'Complete'.",
              " Setting 'type' to 'Complete'")
      type <- "Complete"
    }
    df <- data.frame(Class = factor(seq_len(G)))
  }

  if (!all(sapply(df, is.factor))) {
    stop("Not all columns in the data.frame 'df' are factors")
  }

  stopifnot(G == nrow(df))

  # Make sure the levels of the factors are ordered correctly
  for (i in seq_len(ncol(df))) {
    df[[i]] <- factor(df[[i]], levels = unique(df[[i]]))
  }

  # Construct penalty matrix class
  if (type == "Complete") {

    M <- matrix("fusion", G, G)
    rownames(M) <- colnames(M) <- Reduce(":", df)
    diag(M) <- "ridge"
    return(M)

  } else if (type == "CartesianEqual" || type == "CartesianUnequal") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.cartesianProd, adj.mats)

    if (type == "CartesianEqual") {
      M[M != ""] <- "fusion"
    }
    diag(M) <- "ridge"
    return(M)

  } else if (type == "TensorProd") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.tensorProd, adj.mats)
    diag(M) <- "ridge"
    return(M)

  } else {

    stop("type =", type, "not implemented yet!")

  }
}




##------------------------------------------------------------------------------
##
## To fuse or not to fuse --- Test H0: Omega_1 = ... = Omega_G
##
##------------------------------------------------------------------------------

.scoreStatistic <- function(Plist, Slist, ns) {
  ##############################################################################
  # - Function for computing the score statistic
  # - Plist > A list of precision matrices
  # - Slist > A list of sample covariance matrices
  # - ns    > A vector with the same length as Plist and Slist of sample sizes
  ##############################################################################

  stopifnot(length(Plist) == length(Slist))
  stopifnot(length(ns) == length(Plist))

  U <- 0
  for (g in seq_along(ns)) {
    X <- ns[g]*(Plist[[g]] - Slist[[g]])
    diag(X) <- 0.5*diag(X)
    U <- U + sum(X * (Plist[[g]] %*% X %*% Plist[[g]]))
  }
  return(U)
}



.scambleYlist <- function(Ylist) {
  ##############################################################################
  # - Function for permuting the class labels of Ylist, equivalent to
  #   scrambling/permuting all obervations.
  # - Ylist > A list of observations matrices for each class
  ##############################################################################

  ns <- sapply(Ylist, nrow)
  cl <- factor(rep(names(Ylist), ns))
  Y  <- as.data.frame(do.call(rbind, Ylist))
  out <- split(Y, sample(cl))
  out <- lapply(out, as.matrix)
  return(out)
}



fused.test <- function(Ylist, Tlist, lambda,
                       n.permutations = 100, verbose = FALSE, ...) {
  ##############################################################################
  # - Function for testing the null hypothesis that all population precision
  #   matrices are equal.
  # - Ylist  > A list of G observation matrices for each class.
  # - Tlist  > A list of G p.d. target matrices.
  # - lambda > The non-negative, symmetric G by G penalty matrix
  # - n.permutations > The number of permutation to perform
  # - verbose        > Print out extra progress information
  # - ...            > Arguments passed to ridgeP.fused
  #
  # NOTES:
  # - The test performed is conditional on the supplied penalties and targets.
  # - Returns a object of class "ptest", which has plot, print, and summary
  #   methods.
  ##############################################################################

  stopifnot(length(Ylist) == length(Tlist))
  stopifnot(nrow(lambda) == length(Ylist))
  stopifnot(ncol(lambda) == length(Ylist))

  G <- length(Ylist)
  ns <- sapply(Ylist, nrow)
  n.tot <- sum(ns)
  lambda.null <- G*lambda/sum(ns)

  # Compute observed statistic
  if (verbose) {message("Computing the observed score statistic... ")}
  Slist <- lapply(Ylist, covML)
  Spool.obs <- pooledS(Slist, ns)
  Plist.obs <- list()
  for (i in seq_len(G)) {
    Plist.obs[[i]] <- .armaRidgeP(Spool.obs, target = Tlist[[i]],
                                  lambda = lambda.null[i,i])
  }

  Uobs <- .scoreStatistic(Plist = Plist.obs, Slist = Slist, ns = ns)

  # Approximate null distribution by permutation
  if (verbose) {message("Computing the score statistics under permutation... ")}
  Unull <- numeric()
  for (j in seq_len(n.permutations)) {
    Ylist.tmp <- .scambleYlist(Ylist)  # Permute class labels
    Spool.tmp <- pooledS(lapply(Ylist.tmp, covML), ns)
    Plist.null <- list()
    for (i in seq_len(G)) {
      Plist.null[[i]] <- .armaRidgeP(Spool.tmp,
                                     target = Tlist[[i]],
                                     lambda = lambda.null[i,i])
    }
    Unull[j] <- .scoreStatistic(Plist = Plist.null,
                                Slist = lapply(Ylist.tmp, covML), ns = ns)

    if (verbose && j %% 10 == 0) {
      cat(sprintf("%d of %d done\n", j, n.permutations))
    }
  }

  # Return results
  ans <- list(observed = Uobs, null.dist = Unull)
  class(ans) <- "ptest"
  return(ans)
}



print.ptest <- function(x, digits = 4L, ...) {
  ##############################################################################
  # - Print function for ptest objects
  # - x > A ptest object. Usually created by fused.test()
  ##############################################################################

  x$n.extreme <- sum(x$null.dist >= x$observed)
  x$n.permutations <- length(x$null.dist)
  x$p.val.unbiased <- x$n.extreme/x$n.permutations
  x$p.val.biased <- (x$n.extreme + 1)/(x$n.permutations + 1)
  pval <- format.pval(x$p.val.unbiased, digits = digits,
                      eps = 1/x$n.permutations)

  cat("\nScore-based permutation test\n\n")
  cat("Null hypothesis: Population precision matrices are equal\n")
  cat("Alternative:     Population precision matrices are not equal\n\n")
  cat(sprintf("Observed statistic: U = %0.3f, ", x$observed))
  cat(sprintf("p-value %s\n", ifelse(grepl("<", pval), pval, paste("=", pval))))
  cat("Summary of null distribution obtained by permutation:\n")
  print(summary(x$null.dist, digits = digits, ...))
  return(invisible(x))
}



summary.ptest <- function(object, ...) {
  ##############################################################################
  # - Summary function for ptest objects
  # - x > A ptest object. Usually created by fused.test()
  ##############################################################################

  object <- print.ptest(object, ...)
  cat("\nThe number of extreme observations under the null hypothesis")
  cat(sprintf("\nwas %d out of %d permutations.",
              object$n.extreme, object$n.permutations))

  return(invisible(object))
}



hist.ptest <- function(x, add.extra = TRUE, ...) {
  ##############################################################################
  # - Plot function for ptest objects as a histogram
  # - x          > A ptest object. Usually created by fused.test()
  # - add.extra  > Add the rug of values under the null distribution and
  #                the observed values? Default is TRUE.
  # - ...        > Arguments passed to hist. See ?hist
  ##############################################################################

  hist.args <- list(...)
  if (!hasArg("xlim")) {
    hist.args$xlim <- range(x$null.dist, x$observed)
  }
  if (!hasArg("col")) {
    hist.args$col <- "gray"
  }
  if (!hasArg("main")) {
    hist.args$main <- "Null distribution of U"
  }
  if (!hasArg("xlab")) {
    hist.args$xlab <- "U"
  }
  out <- do.call(hist, c(list(x$null.dist), hist.args), quote = FALSE)
  out$xname <- "x$null.dist"

  if (add.extra) {
    rug(x$null.dist)
    abline(v = x$observed, col = "red", lwd = 2)
    text(x$observed, y = par()$usr[4],
         labels = "Observed U", pos = 3, xpd = TRUE)
  }
  return(invisible(out))
}



plot.ptest <- function(x, add.extra = TRUE, ...) {
  ##############################################################################
  # - Alias for plot.ptest
  ##############################################################################

  hist.ptest(x, add.extra = add.extra, ...)
}




##------------------------------------------------------------------------------
##
## Sparsification and network stats
##
##------------------------------------------------------------------------------

sparsify.fused <- function(Plist, ...) {
  ##############################################################################
  # - Simple wrapper for sparsify. See help(sparsify).
  # - Plist > A list of precision matrices.
  # - ...   > Arguments passed to sparsify.
  ##############################################################################

  return(lapply(Plist, sparsify, ...))
}



GGMnetworkStats.fused <- function(Plist) {
  ##############################################################################
  # - Simple wrapper for GGMnetworkStats. See help(GGMnetworkStats).
  # - Plist > A list of sparse precision matrices.
  ##############################################################################

  res <- lapply(Plist, GGMnetworkStats, as.table = TRUE)
  if (is.null(names(res))) {
    names(res) <- seq_along(Plist)
  }
  return(as.data.frame(res))
}



GGMpathStats.fused <- function(sparsePlist, ...) {
  ##############################################################################
  # - A wrapper for GGMpathStats in the fused case. See GMMpathStats.
  # - sparsePlist > A list of sparsified precision matrices
  # - ...         > Arguments passed to GGMpathStats
  ##############################################################################

  # See if verbose is in ... and set to GGMpathStats default if not
  args <- list(...)
  if (is.null(args[["verbose"]])) {
    verbose <- formals(GGMpathStats)$verbose
  }

  # Run through each class
  res <- vector("list", length(sparsePlist))
  names(res) <- names(sparsePlist)
  for (g in seq_along(res)) {
    if (verbose) {
      cat("\n\n========================================\n",
          "Class: ", names(res)[g], "\n", sep = "")
    }
    res[[g]] <- GGMpathStats(sparsePlist[[g]], ...)
  }
  return(res)
}




