
##' Summary for phylo4/phylo4d objects
##'
##' Summary of information for the tree (\code{phylo4} only) and/or the
##' associated data (\code{phylo4d}).
##'
##' @name summary-methods
##' @docType methods
##' @param object a phylo4d object
##' @param quiet Should the summary be displayed on screen?
##' @param \dots optional additional elements (not in use)
##'
##' @return The \code{nodeType} method returns named vector which has
##' the type of node (internal, tip, root) for value, and the node number
##' for name
##'
##' The \code{summary} method invisibly returns a list with the
##' following components: \item{list("name")}{the name of the object}
##'
##' \item{list("nb.tips")}{the number of tips}
##'
##' \item{list("nb.nodes")}{the number of nodes}
##'
##' \item{list("mean.el")}{mean of edge lengths}
##'
##' \item{list("var.el")}{variance of edge lengths (estimate for population) }
##'
##' \item{list("sumry.el")}{summary (i.e. range and quartiles) of the
##' edge lengths}
##'
##' \item{list("degree")}{(optional) type of polytomy for each node:
##' \sQuote{node}, \sQuote{terminal} (all descendants are tips) or
##' \sQuote{internal} (at least one descendant is an internal node);
##' displayed only when there are polytomies}
##'
##' \item{list("sumry.tips")}{(optional) summary for the data
##' associated with the tips}
##'
##' \item{list("sumry.nodes")}{(optional) summary for the data
##' associated with the internal nodes}
##'
##' @author Ben Bolker, Thibaut Jombart, Francois Michonneau
##' @seealso \code{\link{phylo4d-methods}} constructor and
##' \code{\linkS4class{phylo4d}} class.
##' @keywords methods
##' @aliases summary
##' @include phylo4-methods.R
##' @include phylo4d-methods.R
##' @exportMethod summary
##' @examples
##'   tOwls <- "(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
##'   tree.owls <- ape::read.tree(text=tOwls)
##'   P1 <- as(tree.owls, "phylo4")
##'   P1
##'   summary(P1)
##'   nodeType(P1)
##'
##'   ## summary of a polytomous tree
##'   E <- matrix(c(
##'       8,  9,
##'       9, 10,
##'      10,  1,
##'      10,  2,
##'       9,  3,
##'       9,  4,
##'       8, 11,
##'      11,  5,
##'      11,  6,
##'      11,  7,
##'       0,  8), ncol=2, byrow=TRUE)
##'
##'   P2 <- phylo4(E)
##'   nodeLabels(P2) <- as.character(nodeId(P2, "internal"))
##'   plot(P2, show.node.label=TRUE)
##'   sumryP2 <- summary(P2)
##'   sumryP2
##'
setGeneric("summary")

##' @rdname summary-methods
##' @aliases summary,phylo4-method
setMethod("summary", signature(object="phylo4"),
  function(object, quiet=FALSE) {

      res <- list()

      ## build the result object
      res$name <- deparse(substitute(object, sys.frame(-1)))
      res$nb.tips <- nTips(object)
      res$nb.nodes <- nNodes(object)

      if(hasEdgeLength(object)) {
          edge.length <- edgeLength(object)
          res$mean.el <- mean(edge.length, na.rm=TRUE)
          res$var.el <- stats::var(edge.length, na.rm=TRUE)
          if (isRooted(object) && is.na(edgeLength(object, rootNode(object)))) {
              root.index <- match(edgeId(object, "root"), names(edge.length))
              res$sumry.el <- summary(edge.length[-root.index])
          } else {
              res$sumry.el <- summary(edge.length)
          }
      }

      ## check for polytomies
      if (hasPoly(object)) {
          E <- edges(object)
          temp <- tabulate(E[,1][!is.na(E[, 1])])
          degree <- temp[E[,1][!is.na(E[, 1])]] # contains the degree of the ancestor for all edges
          endsAtATip <- !(E[,2] %in% E[,1])
          terminPoly <- (degree>2) & endsAtATip
          internPoly <- (degree>2) & !endsAtATip
          res$degree <- degree
          res$polytomy <- rep("none",nrow(E))
          res$polytomy[terminPoly] <- "terminal"
          res$polytomy[internPoly] <- "internal"
          ## now just keep information about nodes (not all edges)
          nod <- unique(E[,1])
          idx <- match(nod,E[,1])
          res$degree <- res$degree[idx]
          names(res$degree) <- nodeLabels(object)
          res$polytomy <- res$polytomy[idx]
          names(res$polytomy) <- nodeLabels(object)
      }

      ## model info
      res$loglik <- attr(object, "loglik")
      res$para <- attr(object, "para")
      res$xi <- attr(object, "xi")

      ## if quiet, stop here
      if(quiet) return(invisible(res))

      ## now, print to screen is !quiet
      cat("\n Phylogenetic tree :", res$name, "\n\n")
      cat(" Number of tips    :", res$nb.tips, "\n")
      cat(" Number of nodes   :", res$nb.nodes, "\n")
      ## cat("  ")
      if(hasEdgeLength(object)) {
          cat(" Branch lengths:\n")
          cat("        mean         :", res$mean.el, "\n")
          cat("        variance     :", res$var.el, "\n")
          cat("        distribution :\n")
          print(res$sumry.el)
      }
      else {
          cat(" Branch lengths    : No branch lengths.\n")
      }
      if (hasPoly(object)) {
          cat("\nDegree of the nodes  :\n")
          print(res$degree)
          cat("\n")
          cat("Types of polytomy:\n")
          print(res$polytomy)
          cat("\n")
      }

      if (!is.null(attr(object, "loglik"))) {
          cat("Phylogeny estimated by maximum likelihood.\n")
          cat("  log-likelihood:", attr(object, "loglik"), "\n\n")
          npart <- length(attr(object, "para"))
          for (i in 1:npart) {
              cat("partition ", i, ":\n", sep = "")
              print(attr(object, "para")[[i]])
              if (i == 1)
                  next
              else cat("  contrast parameter (xi):", attr(object,"xi")[i - 1], "\n")
        }
      }
      return(invisible(res))

  })


##' @rdname summary-methods
##' @aliases summary,phylo4d-method
setMethod("summary", signature(object="phylo4d"),
 function(object, quiet=FALSE) {
    x <- object
    res <- summary(as(x, "phylo4"), quiet=quiet)
    res$name <- deparse(substitute(object, sys.frame(-1)))
    tips <- tdata(object, "tip")
    nodes <- tdata(object, "internal")

    if (!quiet)
        cat("\nComparative data:\n")

    if (nrow(tips) > 0) {
        if(!quiet) {
            cat("\nTips: data.frame with", nTips(object), "taxa and",
                ncol(tips), "variable(s) \n\n")
        }
        sumry.tips <- summary(tips)
        res$sumry.tips <- sumry.tips
        if (!quiet)
            print(sumry.tips)
    }
    else {
        if (!quiet)
            cat("\nObject contains no tip data.")
    }
    if (nrow(nodes) > 0) {
        if (!quiet) {
            cat("\nNodes: data.frame with", nNodes(object), "internal nodes and",
                ncol(nodes), "variables \n\n")
        }
        sumry.nodes <- summary(nodes)
        res$sumry.nodes <- sumry.nodes
        if (!quiet)
            print(sumry.nodes)
    }
    else {
        if(!quiet)
            cat("\nObject contains no node data.\n")
    }
    invisible(res)
})

##' @rdname summary-methods
##' @aliases nodeType
##' @export
setGeneric("nodeType", function(object) {
    standardGeneric("nodeType")
})

##' @rdname summary-methods
##' @aliases nodeType,phylo4-method
setMethod("nodeType", signature(object="phylo4"),
 function(object) {
    if(nTips(object) == 0)
        return(NULL)
    else {
        ## strip out the root ancestor
        nodesVect <- as.vector(edges(object))
        nodesVect <- nodesVect[nodesVect != 0]
        ## get a sorted list of the unique nodes
        listNodes <- sort(unique(nodesVect))
        t <- rep("internal", length(listNodes)) # FM: internal is default (I think it's safer)
        names(t) <- listNodes

        ## node number of real internal nodes
        iN <- names(table(edges(object)[,1]))
        ## node number that are not internal nodes (ie that are tips)
        tN <- names(t)[!names(t) %in% iN]
        t[tN] <- "tip"

        ## if the tree is rooted
        if(isRooted(object)) t[rootNode(object)] <- "root"

        return(t)
    }
})
