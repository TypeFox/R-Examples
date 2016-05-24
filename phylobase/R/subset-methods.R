################
## subset phylo4
################

##' Methods for creating subsets of phylogenies
##'
##' Methods for creating subsets of phylogenies, based on pruning a
##' tree to include or exclude a set of terminal taxa, to include all
##' descendants of the MRCA of multiple taxa, or to return a subtree
##' rooted at a given node.
##'
##' The \code{subset} methods must be called using no more than one of
##' the four main subsetting criteria arguments (\code{tips.include},
##' \code{tips.exclude}, \code{mrca}, or \code{node.subtree}).  Each
##' of these arguments can be either character or numeric.  In the
##' first case, they are treated as node labels; in the second case,
##' they are treated as node numbers.  For the first two arguments,
##' any supplied tips not found in the tree (\code{tipLabels(x)}) will
##' be ignored, with a warning.  Similarly, for the \code{mrca}
##' argument, any supplied tips or internal nodes not found in the
##' tree will be ignored, with a warning.  For the \code{node.subtree}
##' argument, failure to provide a single, valid internal node will
##' result in an error.
##'
##' Although \code{prune} is mainly intended as the workhorse function
##' called by \code{subset}, it may also be called directly.  In
##' general it should be equivalent to the \code{tips.exclude} form of
##' \code{subset} (although perhaps with less up-front error
##' checking).
##'
##' The "[" operator, when used as \code{x[i]}, is similar to the
##' \code{tips.include} form of \code{subset}.  However, the indices
##' used with this operator can also be logical, in which case the
##' corresponding tips are assumed to be ordered as in \code{nodeId(x,
##' "tip")}, and recycling rules will apply (just like with a vector
##' or a matrix).  With a \linkS4class{phylo4d} object 'x',
##' \code{x[i,j]} creates a subset of \code{x} taking \code{i} for a
##' tip index and \code{j} for the index of data variables in
##' \code{tdata(geospiza, "all")}.  Note that the second index is
##' optional: \code{x[i, TRUE]}, \code{x[i,]}, and \code{x[i]} are all
##' equivalent.
##'
##' Regardless of which approach to subsetting is used, the argument
##' values must be such that at least two tips are retained.
##'
##' If the most recent common ancestor of the retained tips is not the
##' original root node, then the root node of the subset tree will be
##' a descendant of the original root.  For rooted trees with non-NA
##' root edge length, this has implications for the new root edge
##' length.  In particular, the new length will be the summed edge
##' length from the new root node back to the original root (including
##' the original root edge).  As an alternative, see the examples for
##' a way to determine the length of the edge that was immediately
##' ancestral to the new root node in the original tree.
##'
##' Note that the correspondance between nodes and labels (and data in
##' the case of \linkS4class{phylo4d}) will be retained after all
##' forms of subsetting.  Beware, however, that the node numbers (IDs)
##' will likely be altered to reflect the new tree topology, and
##' therefore cannot be compared directly between the original tree
##' and the subset tree.
##'
##' @name subset-methods
##' @docType methods
##' @param x an object of class \code{"phylo4"} or \code{"phylo4d"}
##' @param tips.include A vector of tips to include in the subset tree
##' @param tips.exclude A vector of tips to exclude from the subset
##' tree
##' @param mrca A vector of nodes for determining the most recent
##' common ancestor, which is then used as the root of the subset tree
##' @param node.subtree A single internal node specifying the root of
##' the subset tree
##' @param trim.internal A logical specifying whether to remove
##' internal nodes that no longer have tip descendants in the subset
##' tree
##' @param i (\code{[} method) An index vector indicating tips to
##' include
##' @param j (\code{[} method, phylo4d only) An index vector
##' indicating columns of node/tip data to include
##' @param drop (not in use: for compatibility with the generic method)
##' @param \dots optional additional parameters (not in use)
##' @return an object of class \code{"phylo4"} or \code{"phylo4d"}
##' @section Methods: \describe{ \item{x = "phylo4"}{subset tree}
##' \item{x = "phylo4d"}{subset tree and corresponding node and tip
##' data} }
##' @author Jim Regetz \email{regetz@@nceas.ucsb.edu}\cr Steven Kembel
##' \email{skembel@@berkeley.edu}\cr Damien de Vienne
##' \email{damien.de-vienne@@u-psud.fr}\cr Thibaut Jombart
##' \email{jombart@@biomserv.univ-lyon1.fr}
##' @keywords methods
##' @exportMethod subset
##' @rdname subset-methods
##' @aliases subset
##' @examples
##' data(geospiza)
##' nodeLabels(geospiza) <- paste("N", nodeId(geospiza, "internal"), sep="")
##' geotree <- extractTree(geospiza)
##'
##' ## "subset" examples
##' tips <- c("difficilis", "fortis", "fuliginosa", "fusca", "olivacea",
##'     "pallida", "parvulus", "scandens")
##' plot(subset(geotree, tips.include=tips))
##' plot(subset(geotree, tips.include=tips, trim.internal=FALSE))
##' plot(subset(geotree, tips.exclude="scandens"))
##' plot(subset(geotree, mrca=c("scandens","fortis","pauper")))
##' plot(subset(geotree, node.subtree=18))
##'
##' ## "prune" examples (equivalent to subset using tips.exclude)
##' plot(prune(geotree, tips))
##'
##' ## "[" examples (equivalent to subset using tips.include)
##' plot(geotree[c(1:6,14)])
##' plot(geospiza[c(1:6,14)])
##'
##' ## for phylo4d, subset both tips and data columns
##' geospiza[c(1:6,14), c("wingL", "beakD")]
##'
##' ## note handling of root edge length:
##' edgeLength(geotree)['0-15'] <- 0.1
##' geotree2 <- geotree[1:2]
##' ## in subset tree, edge of new root extends back to the original root
##' edgeLength(geotree2)['0-3']
##' ## edge length immediately ancestral to this node in the original tree
##' edgeLength(geotree, MRCA(geotree, tipLabels(geotree2)))
setGeneric("subset")

##' @rdname subset-methods
##' @aliases subset,phylo4-method
setMethod("subset", "phylo4", function(x, tips.include=NULL,
    tips.exclude=NULL, mrca=NULL, node.subtree=NULL, ...) {
    ##  FIXME: could eliminate NULL and make the test
    ## if (!missing) rather than if (!is.null)
    ## (might have to change next line?)
    if (sum(!sapply(list(tips.include, tips.exclude, mrca,
        node.subtree), is.null))>1) {
        stop("must specify at most one criterion for subsetting")
    }
    all.tips <- nodeId(x, "tip")
    if (!is.null(tips.include)) {
        nodes <- getNode(x, tips.include, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        kept <- nodes[is.valid.tip]
        dropped <- setdiff(all.tips, kept)
        unknown <- tips.include[!is.valid.tip]
    } else if (!is.null(tips.exclude)) {
        nodes <- getNode(x, tips.exclude, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        dropped <- nodes[is.valid.tip]
        kept <- setdiff(all.tips, dropped)
        unknown <- tips.exclude[!is.valid.tip]
    } else if (!is.null(mrca)) {
        nodes <- getNode(x, mrca, missing="OK")
        is.valid.node <- nodes %in% nodeId(x, "all")
        mnode <- MRCA(x, nodes[is.valid.node])
        if (length(mnode)!=1) {
            stop("mrca must include at least one valid node")
        }
        kept <- descendants(x, mnode)
        dropped <- setdiff(all.tips, kept)
        unknown <- mrca[!is.valid.node]
    } else if (!is.null(node.subtree)) {
        node <- getNode(x, node.subtree, missing="OK")
        if (length(node)!=1 || !(node %in% nodeId(x, "internal"))) {
            stop("node.subtree must be a single valid internal node")
        }
        kept <- descendants(x, node)
        dropped <- setdiff(all.tips, kept)
        unknown <- numeric(0)
    } else {
        kept <- getNode(x, nodeId(x, "tip"))
        dropped <- numeric(0)
        unknown <- numeric(0)
    }
    if (length(unknown)>0) {
        warning("invalid nodes ignored: ", paste(unknown,
            collapse=", "))
    }
    if (length(kept)<2) {
        stop("0 or 1 tips would remain after subsetting")
    }
    if (length(dropped)==0) return(x)
    return(prune(x, dropped, ...))
})

###############
# '[' operator
###############

## Consider using some combination of these for stricter argument
## checking? Not implementing now because extra arguments are just
## ignored, which is fairly common S4 method behavior:
## * in "[" methods for phylo4:
##    if (nargs()>2) stop("unused arguments")
## * in "[" methods for both phylo4 and phylo4d:
##    if (!missing(...)) stop("unused argument(s)")

##' @rdname subset-methods
##' @aliases [
##' @exportMethod [
setGeneric("[")

##### --------  phylo4 '[' methods

##' @rdname subset-methods
##' @aliases [,phylo4,character,missing-method
setMethod("[", signature(x="phylo4", i="character", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})

##' @rdname subset-methods
##' @aliases [,phylo4,numeric,missing-method
setMethod("[", signature(x="phylo4", i="numeric", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})

##' @rdname subset-methods
##' @aliases [,phylo4,logical,missing-method
setMethod("[", signature(x="phylo4", i="logical", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=nodeId(x, "tip")[i])
})

##' @rdname subset-methods
##' @aliases [,phylo4,missing,missing-method
setMethod("[", signature(x="phylo4", i="missing", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    return(x)
})

##### --------  phylo4d '[' methods

##' @rdname subset-methods
##' @aliases [,phylo4d,ANY,character,missing-method
setMethod("[", signature(x="phylo4d", i="ANY", j="character",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})

##' @rdname subset-methods
##' @aliases [,phylo4d,ANY,numeric,missing-method
setMethod("[", signature(x="phylo4d", i="ANY", j="numeric",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})

##' @rdname subset-methods
##' @aliases [,phylo4d,ANY,logical,missing-method
setMethod("[", signature(x="phylo4d", i="ANY", j="logical",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})

## borrow from Matrix package approach of trapping invalid usage
##' @rdname subset-methods
##' @aliases [,phylo4,ANY,ANY,ANY-method
setMethod("[", signature(x="phylo4", i="ANY", j="ANY", drop="ANY"),
    function(x, i, j, ..., drop) {
    stop("invalid argument(s)")
})

##### -------- prune

##' @rdname subset-methods
##' @aliases prune
##' @export
setGeneric("prune", function(x, ...) {
    standardGeneric("prune")
})

## return characters, sorted in NUMERIC order
.chnumsort <- function(x) {
  as.character(sort(as.numeric(x)))
}

##' @rdname subset-methods
##' @aliases prune,phylo4-method
setMethod("prune", "phylo4",
    function(x, tips.exclude, trim.internal=TRUE) {

    makeEdgeNames <- function(edge) {
        paste(edge[,1], edge[,2], sep="-")
    }

    ## drop tips and obsolete internal nodes from edge matrix
    tip.drop <- getNode(x, tips.exclude, missing="fail")
    tip.keep <- setdiff(nodeId(x, "tip"), tip.drop)
    nodes <- nodeId(x, "all")
    node.keep <- rep(FALSE, length(nodes))
    node.keep[tip.keep] <- TRUE
    if (trim.internal) {
        if (edgeOrder(x) == "postorder") {
            edge.post <- edges(x)
        } else {
            edge.post <- edges(reorder(x, "postorder"))
        }
        for (i in seq_along(edge.post[,2])) {
            if (node.keep[edge.post[i,2]]) {
                node.keep[edge.post[i,1]] <- TRUE
            }
        }
    } else {
        node.keep[nodeId(x, "internal")] <- TRUE
    }
    edge.new <- edges(x)[edges(x)[,2] %in% nodes[node.keep], ]

    ## remove singletons
    edge.length.new <- edgeLength(x)
    edge.label.new <- edgeLabels(x)
    singletons <- which(tabulate(edge.new[edge.new[, 1] != 0, 1])==1)
    while (length(singletons)>0) {
        sing.node <- singletons[1]

        ## update edge matrix
        edges.drop <- which(edge.new==sing.node, arr.ind=TRUE)[,"row"]
        sing.edges <- edge.new[edges.drop,]
        edge.new[edges.drop[2], ] <- c(sing.edges[2,1], sing.edges[1,2])
        edge.new <- edge.new[-edges.drop[1], ]

        ## update edge lengths and edge labels
        edge.names.drop <- makeEdgeNames(sing.edges)
        edge.name.new <- paste(sing.edges[2,1], sing.edges[1,2], sep="-")
        edge.length.new[edge.name.new] <-
            sum(edge.length.new[edge.names.drop])
        edge.length.new <- edge.length.new[-match(edge.names.drop,
            names(edge.length.new))]
        edge.label.new[edge.name.new] <- NA
        edge.label.new <- edge.label.new[-match(edge.names.drop,
            names(edge.label.new))]

        singletons <- which(tabulate(edge.new[edge.new[, 1] != 0, 1])==1)
    }

    ## remove dropped elements from labels
    label.new <- labels(x)[names(labels(x)) %in% edge.new]

    ## subset and order edge.length and edge.label with respect to edge
    edge.names <- makeEdgeNames(edge.new)
    edge.length.new <- edge.length.new[edge.names]
    edge.label.new <- edge.label.new[edge.names]

    if (!trim.internal) {
        ## make sure now-terminal internal nodes are treated as tips
        tip.now <- setdiff(edge.new[,2], edge.new[,1])
        tip.add <- tip.now[tip.now>nTips(x)]
        if (length(tip.add)>0) {
            ind <- match(tip.add, names(label.new))

            ## node renumbering workaround to satisfy plot method
            newid <- sapply(tip.add, function(tip) descendants(x, tip)[1])
            names(label.new)[ind] <- newid
            edge.new[match(tip.add, edge.new)] <- newid
            tip.now[match(tip.add, tip.now)] <- newid

            isTip <- edge.new %in% tip.now
            edge.new[isTip] <- match(edge.new[isTip],
            sort(unique.default(edge.new[isTip])))
        }
    }

    ## renumber nodes in the edge matrix
    edge.new[] <- match(edge.new, sort(unique.default(edge.new))) - 1L

    ## update corresponding element names in the other slots
    edge.names <- makeEdgeNames(edge.new)
    names(edge.length.new) <- edge.names
    names(edge.label.new) <- edge.names
    label.new <- label.new[order(as.numeric(names(label.new)))]
    names(label.new) <- seq_along(label.new)

    ## update, check, then return the pruned phylo4 object
    x@edge <- edge.new
    ##TODO would prefer to leave out NA from edge.length slot, but can't
    x@edge.length <- edge.length.new
    x@edge.label <- edge.label.new[!is.na(edge.label.new)]
    x@label <- label.new[!is.na(label.new)]
    if(is.character(checkval <- checkPhylo4(x))) {
        stop(checkval)
    } else {
        return(x)
    }

})

##' @rdname subset-methods
##' @aliases prune,phylo4d-method
setMethod("prune", "phylo4d",
    function(x, tips.exclude, trim.internal=TRUE) {

    tree <- extractTree(x)
    phytr <- prune(tree, tips.exclude, trim.internal)

    ## create temporary phylo4 object with complete and unique labels
    tmpLbl <- .genlab("n", nTips(x)+nNodes(x))
    tmpPhy <- tree
    labels(tmpPhy, "all") <- tmpLbl
    tmpPhytr <- prune(tmpPhy, getNode(x, tips.exclude), trim.internal)

    ## get node numbers to keep
    oldLbl <- labels(tmpPhy, "all")
    newLbl <- labels(tmpPhytr, "all")
    wasKept <- oldLbl %in% newLbl
    nodesToKeep <- as.numeric(names(oldLbl[wasKept]))

    ## subset original data, and update names
    allDt <- x@data[match(nodesToKeep, rownames(x@data)), , drop=FALSE]
    rownames(allDt) <- match(newLbl, oldLbl[wasKept])

    phytr <- phylo4d(phytr, all.data=allDt, match.data=TRUE)

    phytr
})

## setMethod("prune","ANY",
##           function(phy, tip, trim.internal = TRUE, subtree = FALSE,
##                    ,...) {
##             if (class(phy)=="phylo") {
##               ape::prune(phy, tip, trim.internal, subtree)
##               } else stop("no prune method available for",
##                     deparse(substitute(phy)),
##                     "(class",class(phy),")")
##           })
