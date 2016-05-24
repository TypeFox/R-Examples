
##' Converting between phylo4/phylo4d and other phylogenetic tree
##' formats
##'
##' Translation functions to convert between phylobase objects
##' (\code{phylo4} or \code{phylo4d}), and objects used by other
##' comparative methods packages in R: \code{ape} objects
##' (\code{phylo}, \code{multiPhylo}), \code{RNeXML} object
##' (\code{nexml}), \code{ade4} objects (\code{phylog}, \emph{now
##' deprecated}), and to \code{data.frame} representation.
##'
##' @name setAs
##' @docType methods
##' @section Usage: \code{as(object, class)}
##' @author Ben Bolker, Thibaut Jombart, Marguerite Butler, Steve
##' Kembel, Francois Michonneau
##' @seealso generic \code{\link[methods]{as}},
##' \code{\link{phylo4-methods}}, \code{\link{phylo4d-methods}},
##' \code{\link{extractTree}}, \code{nexml} class from the
##' \code{RNeXML} package, \code{\link[ade4]{phylog}} from the
##' \code{ade4} package and \code{\link[ape]{as.phylo}} from the
##' \code{ape} package.
##' @keywords methods
##' @rdname setAs-methods
##' @aliases as as-method as,phylo,phylo4-method
##' @include phylo4-methods.R
##' @include phylo4d-methods.R
##' @include oldclasses-class.R
##' @examples
##' tree_string <- "(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
##' tree.owls <- ape::read.tree(text=tree_string)
##' ## round trip conversion
##' tree_in_phylo <- tree.owls                  # tree is a phylo object
##' (tree_in_phylo4 <- as(tree.owls,"phylo4"))  # phylo converted to phylo4
##' identical(tree_in_phylo,as(tree_in_phylo4,"phylo"))
##' ## test if phylo, and phylo4 converted to phylo are identical
##' ## (no, because of dimnames)
##'
##' ## Conversion to phylog (ade4)
##' as(tree_in_phylo4, "phylog")
##'
##' ## Conversion to data.frame
##' as(tree_in_phylo4, "data.frame")
##'
##' ## Conversion to phylo (ape)
##' as(tree_in_phylo4, "phylo")
##'
##' ## Conversion to phylo4d, (data slots empty)
##' as(tree_in_phylo4, "phylo4d")
setAs("phylo", "phylo4", function(from, to) {
  ## fixme SWK kludgy fix may not work well with unrooted trees
  ## TODO should we also attempt to get order information?
  ## BMB horrible kludge to avoid requiring ape explicitly
    ape_is.rooted <- function(phy) {
        if (!is.null(phy$root.edge))
            TRUE
        else if (tabulate(phy$edge[, 1])[length(phy$tip.label) + 1] > 2)
            FALSE
        else TRUE
    }
    if (ape_is.rooted(from)) {
        tip.idx <- 1:nTips(from)
        if (nTips(from) < nrow(from$edge)) {
            int.idx <- (nTips(from)+1):dim(from$edge)[1]
        } else {
            int.idx <- NULL
        }
        root.node <- as.numeric(setdiff(unique(from$edge[,1]),
                                        unique(from$edge[,2])))

        from$edge <- rbind(from$edge[tip.idx,],c(0,root.node),
                           from$edge[int.idx,])
        if (!is.null(from$edge.length)) {
            if (is.null(from$root.edge)) {
                from$edge.length <- c(from$edge.length[tip.idx],
                                      as.numeric(NA),from$edge.length[int.idx])
            }
            else {
                from$edge.length <- c(from$edge.length[tip.idx],
                                      from$root.edge,from$edge.length[int.idx])
            }
        }
        if (!is.null(from$edge.label)) {
            from$edge.label <- c(from$edge.label[tip.idx], NA,
                                 from$edge.label[int.idx])
        }
    }
    newobj <- phylo4(from$edge, from$edge.length, unname(from$tip.label),
                     node.label = from$node.label,
                     edge.label = from$edge.label,
                     order = "unknown")
    oldorder <- attr(from,"order")
    neworder <- if (is.null(oldorder)) {
        "unknown"
    } else if (!oldorder %in% phylo4_orderings) {
        stop("unknown ordering '", oldorder, "' in ape object")
    } else if (oldorder == "cladewise" || oldorder == "preorder") {
        "preorder"
    } else if (oldorder == "pruningwise" || oldorder == "postorder") {
        "postorder"
    }
    if (isRooted(newobj)) {
        if (neworder == "preorder") {
            newobj <- reorder(newobj, order="preorder")
        }
        if (neworder == "postorder") {
            newobj <- reorder(newobj, order="postorder")
        }
    }
    newobj@order <- neworder

    attr(from,"order") <- NULL

    attribs <- attributes(from)
    attribs$names <- NULL
    knownattr <- c("logLik", "origin", "para", "xi")
    known <- names(attribs)[names(attribs) %in% knownattr]
    unknown <- names(attribs)[!names(attribs) %in% c(knownattr, "class", "names")]
    if (length(unknown) > 0) {
        warning(paste("unknown attributes ignored: ", unknown, collapse = " "))
    }
    for (i in known) attr(newobj, i) <- attr(from, i)
    newobj
})

##' @name setAs
##' @rdname setAs-methods
##' @aliases as,phylo,phylo4d-method
setAs("phylo", "phylo4d", function(from, to) {
    phylo4d(as(from, "phylo4"), tip.data = data.frame())
})

##' @name setAs
##' @rdname setAs-methods
##' @aliases as,nexml,phylo4-method
setAs("nexml", "phylo4", function(from, to) {
    phylo4(from)
})

##' @name setAs
##' @rdname setAs-methods
##' @aliases as,nexml,phylo4d-method
setAs("nexml", "phylo4d", function(from, to) {
    phylo4d(from)
})

#######################################################
## Exporting to ape


## BMB: adding an explicit as method, and the warning,
##  here is a very bad idea, because
##   even implicit conversions from phylo4d to phylo4 (e.g.
##  to use inherited methods) will produce the warning

## setAs("phylo4d", "phylo4",function(from,to) {
##   warning("losing data while coercing phylo4d to phylo")
##   phylo4(from@edge, from@edge.length, from@tip.label,
##         from@node.label,from@edge.label,from@order)
## })

##' @name setAs
##' @rdname setAs-methods
##' @aliases as,phylo4,phylo-method
setAs("phylo4", "phylo", function(from, to) {

    if(is.character(checkval <- checkPhylo4(from))) {
      stop(checkval)
    }

    if (inherits(from, "phylo4d"))
        warning("losing data while coercing phylo4d to phylo")

    phy <- list()

    ## Edge matrix (dropping root edge if it exists)
    edgemat <- unname(edges(from, drop.root=TRUE))
    storage.mode(edgemat) <- "integer"
    phy$edge <- edgemat

    ## Edge lengths
    if(hasEdgeLength(from)) {
        edge.length <- edgeLength(from)
        if(isRooted(from)) {
            iRoot <- match(edgeId(from, "root"), names(edge.length))
            phy$edge.length <- unname(edge.length[-iRoot])
        }
        else {
            phy$edge.length <- unname(edge.length)
        }
    }

    ## Tip labels
    phy$tip.label <- unname(tipLabels(from))

    ## nNodes
    phy$Nnode <- as.integer(nNodes(from))

    ## Node labels
    if(hasNodeLabels(from)) {
        phy$node.label <- unname(nodeLabels(from))
    }

    ## Root edge
    if(isRooted(from) && hasEdgeLength(from)) {
        root.edge <- unname(edgeLength(from,rootNode(from)))
        if(!is.na(root.edge)) {
            phy$root.edge <- root.edge
        }
    }

    ## Converting to class phylo
    class(phy) <- "phylo"

    ## Tree order
    ## TODO postorder != pruningwise -- though quite similar
    if (edgeOrder(from) == "unknown") {
        warning("trees with unknown order may be",
                " unsafe in ape")
    }
    else {
        attr(phy, "order") <- switch(edgeOrder(from),
                                     postorder = "unknown",
                                     preorder = "cladewise",
                                     pruningwise = "pruningwise")
    }
    phy
})


## BMB: redundant????
## JR: updated (but untested) to reflect slot changes, in case this ever
##     needs to come out of its commented hibernation
## setAs("phylo4d", "phylo", function(from, to) {
##     y <- list(edge = edges(from, drop.root=TRUE),
##         Nnode = nNodes(from), tip.label = tipLabels(from))
##     class(y) <- "phylo"
##     if (hasEdgeLength(from))
##         y$edge.length <- edgeLength(from)
##     if (hasNodeLabels(from))
##         y$node.label <- nodeLabels(from)
##     #if (!is.na(from@root.edge))
##     #    y$root.edge <- from@root.edge
##    warning("losing data while coercing phylo4d to phylo")
##    y
##})


#######################################################
## Exporting to ade4

##' @name setAs
##' @rdname setAs-methods
##' @aliases setAs,phylo4,phylog-method
setAs("phylo4", "phylog", function(from, to) {
    x <- as(from, "phylo")
    xstring <- write.tree(x, file = "")
    warning("ade4::phylog objects are deprecated, please use the adephylo package instead")
    ade4::newick2phylog(xstring)
})

#######################################################
## Exporting to dataframe

.phylo4ToDataFrame <- function(from, edgeOrder=c("pretty", "real")) {

    edgeOrder <- match.arg(edgeOrder)

    ## Check the phylo4
    if (is.character(checkval <- checkPhylo4(from)))
        stop(checkval)

    ## The order of 'node' defines the order of all other elements
    if (edgeOrder == "pretty") {
        node <- nodeId(from, "all")
        ancestr <- ancestor(from, node)

        # ancestor returns an NA, replace this w/ 0 to construct names correctly
        ancestr[is.na(ancestr)] <- as.integer(0)
    } else {
        E <- edges(from)
        node <- E[, 2]
        ancestr <- E[, 1]
    }

    ## extract and reorder (as needed) other object slots
    nmE <- paste(ancestr, node, sep="-")
    edge.length <- edgeLength(from)
    edge.length <- edge.length[match(nmE, names(edge.length))]

    ndType <- nodeType(from)
    ndType <- ndType[match(node, names(ndType))]
    label <- labels(from, type="all")
    label <- label[match(node, names(label))]

    tDf <- data.frame(label, node, ancestor=ancestr, edge.length,
                    node.type=ndType, row.names=node)
    tDf$label <- as.character(tDf$label)

    if (class(from) == "phylo4d") {
        dat <- tdata(from, "all", label.type="column") # get data

        ## reorder data to edge matrix order, drop labels (first column)
        if(nrow(dat) > 0 && ncol(dat) > 1) {
            dat <- dat[match(rownames(tDf), rownames(dat)), ]
            tDf <- cbind(tDf, dat[ ,-1 , drop=FALSE])
        }
        else {
            cat("No data associated with the tree\n")
       }
    }
    tDf
}

##' @name setAs
##' @rdname setAs-methods
##' @aliases setAs,phylo4,data.frame-method
setAs(from = "phylo4", to = "data.frame", def=function(from) {
    d <- .phylo4ToDataFrame(from, edgeOrder="pretty")
    d
})
