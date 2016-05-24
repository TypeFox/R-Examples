######################
## phylo4d constructor
######################

## TEST ME
## '...' recognized args for data are tipdata and nodedata.
## other recognized options are those known by the phylo4 constructor


##' Combine a phylogenetic tree with data
##'
##' \code{phylo4d} is a generic constructor which merges a
##' phylogenetic tree with data frames to create a combined object of
##' class \code{phylo4d}
##'
##' You can provide several data frames to define traits associated
##' with tip and/or internal nodes. By default, data row names are
##' used to link data to nodes in the tree, with any number-like names
##' (e.g., \dQuote{10}) matched against node ID numbers, and any
##' non-number-like names (e.g., \dQuote{n10}) matched against node
##' labels. Alternative matching rules can be specified by passing
##' additional arguments (listed in the Details section); these
##' include positional matching, matching exclusively on node labels,
##' and matching based on a column of data rather than on row
##' names.
##'
##' Matching rules will apply the same way to all supplied data
##' frames.  This means that you need to be consistent with the row
##' names of your data frames.  It is good practice to use tip and
##' node labels (or node numbers if you use duplicated labels) when
##' you combine data with a tree.
##'
##' If you provide both \code{tip.data} and \code{node.data}, the
##' treatment of columns with common names will depend on the
##' \code{merge.data} argument. If TRUE, columns with the same name in
##' both data frames will be merged; when merging columns of different
##' data types, coercion to a common type will follow standard R
##' rules. If \code{merge.data} is FALSE, columns with common names
##' will be preserved independently, with \dQuote{.tip} and
##' \dQuote{.node} appended to the names. This argument has no effect
##' if \code{tip.data} and \code{node.data} have no column names in
##' common.
##'
##' If you provide \code{all.data} along with either of
##' \code{tip.data} and \code{node.data}, it must have distinct column
##' names, otherwise an error will result. Additionally, although
##' supplying columns with the same names \emph{within} data frames is
##' not illegal, automatic renaming for uniqeness may lead to
##' surprising results, so this practice should be avoided.
##'
##' @name phylo4d-methods
##' @aliases phylo4d
##' @param x an object of class \code{phylo4}, \code{phylo},
##' \code{nexml} or a matrix of edges (see above)
##' @param tip.data a data frame (or object to be coerced to one)
##' containing only tip data (Optional)
##' @param node.data a data frame (or object to be coerced to one)
##' containing only node data (Optional)
##' @param all.data a data frame (or object to be coerced to one)
##' containing both tip and node data (Optional)
##' @param merge.data if both \code{tip.data} and \code{node.data} are
##' provided, should columns with common names will be merged together
##' (default TRUE) or not (FALSE)? See details.
##' @param metadata any additional metadata to be passed to the new object
##' @param edge.length Edge (branch) length. (Optional)
##' @param tip.label A character vector of species names (names of
##' "tip" nodes).  (Optional)
##' @param node.label A character vector of internal node
##' names. (Optional)
##' @param edge.label A character vector of edge (branch)
##' names. (Optional)
##' @param order character: tree ordering (allowable values are listed
##' in \code{phylo4_orderings}, currently "unknown", "preorder"
##' (="cladewise" in \code{ape}), and "postorder", with "cladewise"
##' and "pruningwise" also allowed for compatibility with \code{ape})
##' @param annote any additional annotation data to be passed to the
##' new object
##' @param check.node.labels if \code{x} is of class \code{phylo}, use
##' either \dQuote{keep} (the default) to retain internal node labels,
##' \dQuote{drop} to drop them, or \dQuote{asdata} to convert them to
##' numeric tree data. This argument is useful if the \code{phylo}
##' object has non-unique node labels or node labels with informative
##' data (e.g., posterior probabilities).
##' @param \dots further arguments to control the behavior of the
##' constructor in the case of missing/extra data and where to look
##' for labels in the case of non-unique labels that cannot be stored
##' as row names in a data frame (see Details).
##' @details This is the list of additional arguments that can be used
##' to control matching between the tree and the data:
##'
##' \itemize{
##'
##' \item{match.data}{(logical) should the rownames of the data frame
##' be used to be matched against tip and internal node identifiers?}
##'
##' \item{rownamesAsLabels}{(logical), should the row names of the
##' data provided be matched only to labels (TRUE), or should any
##' number-like row names be matched to node numbers (FALSE and
##' default)}
##'
##' \item{label.type}{character, \code{rownames} or \code{column}:
##' should the labels be taken from the row names of \code{dt} or from
##' the \code{label.column} column of \code{dt}?}
##'
##' \item{label.column}{iff \code{label.type=="column"}, column
##' specifier (number or name) of the column containing tip labels}
##'
##' \item{missing.data}{action to take if there are missing data or if
##' there are data labels that don't match}
##'
##' \item{extra.data}{action to take if there are extra data or if
##' there are labels that don't match}
##'
##' \item{keep.all}{(logical), should the returned data have rows for
##' all nodes (with NA values for internal rows when type='tip', and
##' vice versa) (TRUE and default) or only rows corresponding to the
##' type argument}
##'
##' }
##'
##' Rules for matching rows of data to tree nodes are determined
##' jointly by the \code{match.data} and \code{rownamesAsLabels}
##' arguments. If \code{match.data} is TRUE, data frame rows will be
##' matched exclusively against tip and node labels if
##' \code{rownamesAsLabels} is also TRUE, whereas any all-digit row
##' names will be matched against tip and node numbers if
##' \code{rownamesAsLabels} is FALSE (the default). If
##' \code{match.data} is FALSE, \code{rownamesAsLabels} has no effect,
##' and row matching is purely positional with respect to the order
##' returned by \code{nodeId(phy, type)}.
##'
##' @return An object of class \linkS4class{phylo4d}.
##' @note Checking on matches between the tree and the data will be
##' done by the validity checker (label matches between data and tree
##' tips, number of rows of data vs. number of nodes/tips/etc.)
##' @section Methods: \describe{ \item{x = "phylo4"}{merges a tree of
##' class \code{phylo4} with a data.frame into a \code{phylo4d}
##' object} \item{x = "matrix"}{merges a matrix of tree edges similar
##' to the edge slot of a \code{phylo4} object (or to \$edge of a
##' \code{phylo} object) with a data.frame into a \code{phylo4d}
##' object} \item{x = "phylo"}{merges a tree of class \code{phylo}
##' with a data.frame into a \code{phylo4d} object } }
##' @author Ben Bolker, Thibaut Jombart, Steve Kembel, Francois
##' Michonneau, Jim Regetz
##' @seealso \code{\link{coerce-methods}} for translation
##' functions. The \linkS4class{phylo4d} class; \linkS4class{phylo4}
##' class and \link{phylo4} constructor.
##' @keywords misc
##' @export
##' @docType methods
##' @rdname phylo4d-methods
##' @include phylo4d-class.R
##' @include oldclasses-class.R
##' @examples
##'
##' treeOwls <- "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);"
##' tree.owls.bis <- ape::read.tree(text=treeOwls)
##' try(phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3)), silent=TRUE)
##' obj <- phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3), match.data=FALSE)
##' obj
##' print(obj)
##'
##' ####
##'
##' data(geospiza_raw)
##' geoTree <- geospiza_raw$tree
##' geoData <- geospiza_raw$data
##'
##' ## fix differences in tip names between the tree and the data
##' geoData <- rbind(geoData, array(, dim = c(1,ncol(geoData)),
##'                   dimnames = list("olivacea", colnames(geoData))))
##'
##' ### Example using a tree of class 'phylo'
##' exGeo1 <- phylo4d(geoTree, tip.data = geoData)
##'
##' ### Example using a tree of class 'phylo4'
##' geoTree <- as(geoTree, "phylo4")
##'
##' ## some random node data
##' rNodeData <- data.frame(randomTrait = rnorm(nNodes(geoTree)),
##'                         row.names = nodeId(geoTree, "internal"))
##'
##' exGeo2 <- phylo4d(geoTree, tip.data = geoData, node.data = rNodeData)
##'
##' ### Example using 'merge.data'
##' data(geospiza)
##' trGeo <- extractTree(geospiza)
##' tDt <- data.frame(a=rnorm(nTips(trGeo)), row.names=nodeId(trGeo, "tip"))
##' nDt <- data.frame(a=rnorm(nNodes(trGeo)), row.names=nodeId(trGeo, "internal"))
##'
##' (matchData1 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=FALSE))
##' (matchData2 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=TRUE))
##'
##' ## Example with 'all.data'
##' nodeLabels(geoTree) <- as.character(nodeId(geoTree, "internal"))
##' rAllData <- data.frame(randomTrait = rnorm(nTips(geoTree) + nNodes(geoTree)),
##' row.names = labels(geoTree, 'all'))
##'
##' exGeo5 <- phylo4d(geoTree, all.data = rAllData)
##'
##' ## Examples using 'rownamesAsLabels' and comparing with match.data=FALSE
##' tDt <- data.frame(x=letters[1:nTips(trGeo)],
##'                   row.names=sample(nodeId(trGeo, "tip")))
##' tipLabels(trGeo) <- as.character(sample(1:nTips(trGeo)))
##' (exGeo6 <- phylo4d(trGeo, tip.data=tDt, rownamesAsLabels=TRUE))
##' (exGeo7 <- phylo4d(trGeo, tip.data=tDt, rownamesAsLabels=FALSE))
##' (exGeo8 <- phylo4d(trGeo, tip.data=tDt, match.data=FALSE))
##'
##' ## generate a tree and some data
##' set.seed(1)
##' p3 <- ape::rcoal(5)
##' dat <- data.frame(a = rnorm(5), b = rnorm(5), row.names = p3$tip.label)
##' dat.defaultnames <- dat
##' row.names(dat.defaultnames) <- NULL
##' dat.superset <- rbind(dat, rnorm(2))
##' dat.subset <- dat[-1, ]
##'
##' ## create a phylo4 object from a phylo object
##' p4 <- as(p3, "phylo4")
##'
##' ## create phylo4d objects with tip data
##' p4d <- phylo4d(p4, dat)
##' ###checkData(p4d)
##' p4d.sorted <- phylo4d(p4, dat[5:1, ])
##' try(p4d.nonames <- phylo4d(p4, dat.defaultnames))
##' p4d.nonames <- phylo4d(p4, dat.defaultnames, match.data=FALSE)
##'
##' \dontrun{
##' p4d.subset <- phylo4d(p4, dat.subset)
##' p4d.subset <- phylo4d(p4, dat.subset)
##' try(p4d.superset <- phylo4d(p4, dat.superset))
##' p4d.superset <- phylo4d(p4, dat.superset)
##' }
##'
##' ## create phylo4d objects with node data
##' nod.dat <- data.frame(a = rnorm(4), b = rnorm(4))
##' p4d.nod <- phylo4d(p4, node.data = nod.dat, match.data=FALSE)
##'
##'
##' ## create phylo4 objects with node and tip data
##' p4d.all1 <- phylo4d(p4, node.data = nod.dat, tip.data = dat, match.data=FALSE)
##' nodeLabels(p4) <- as.character(nodeId(p4, "internal"))
##' p4d.all2 <- phylo4d(p4, all.data = rbind(dat, nod.dat), match.data=FALSE)
setGeneric("phylo4d", function(x, ...) { standardGeneric("phylo4d")} )

## first arg is a phylo4
##' @rdname phylo4d-methods
##' @aliases phylo4d,phylo4,phylo4-method
setMethod("phylo4d", "phylo4",
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   merge.data=TRUE, metadata = list(), ...) {
    ## coerce tree to phylo4d
    res <- as(x, "phylo4d")

    ## apply formatData to ensure data have node number rownames and
    ## correct dimensions
    tip.data <- formatData(phy=x, dt=tip.data, type="tip", ...)
    node.data <- formatData(phy=x, dt=node.data, type="internal", ...)
    all.data <- formatData(phy=x, dt=all.data, type="all", ...)

    ## add any data
    res@data <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
        all.data=all.data, merge.data=merge.data)
    ## add any metadata
    res@metadata <- metadata
    return(res)
})


### first arg is a matrix of edges
##' @rdname phylo4d-methods
##' @aliases phylo4d,matrix,matrix-method
setMethod("phylo4d", "matrix",
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   merge.data=TRUE, metadata=list(), edge.length=NULL,
                   tip.label=NULL, node.label=NULL, edge.label=NULL,
                   order="unknown", annote=list(), ...) {
    tree <- phylo4(x, edge.length=edge.length, tip.label=tip.label,
        node.label=node.label, edge.label=edge.label, order=order,
        annote=annote)
    res <- phylo4d(tree, tip.data, node.data, all.data,
        merge.data=merge.data, metadata=metadata, ...)
    return(res)
})


label_to_data <- function(nlab.data, ...) {
    ## convert number-like labels to numeric, other keep as it is
    nlab.data.test <- gsub("[0-9]|\\.", "", nlab.data[!is.na(nlab.data)])
    if (all(nchar(nlab.data.test) == 0 )) {
        nlab.data <- data.frame(labelValues=as.numeric(nlab.data), ...)
    }
    else {
        nlab.data <- data.frame(labelValues=nlab.data, ...)
    }
    nlab.data
}

### first arg is a phylo
##' @rdname phylo4d-methods
##' @aliases phylo4d,phylo,phylo-method
setMethod("phylo4d", "phylo",
          function(x, tip.data=NULL,
                   node.data=NULL, all.data=NULL,
                   check.node.labels=c("keep", "drop", "asdata"),
                   annote=list(), metadata=list(), ...) {

    check.node.labels <- match.arg(check.node.labels)

    if (check.node.labels == "asdata") {
        # FIXME? use.node.names=TRUE won't work with this option b/c
        # node labels are dropped; assumes node.data (if any), phylo
        # node.label, and phylo4 internal nodes are in the same order?

        nlab.data <- x$node.label
        x$node.label <- NULL
        nlab.data[!nzchar(nlab.data)] <- NA

        ## convert number-like labels to numeric, other keep as it is
        nlab.data <- label_to_data(nlab.data)

        tree <- phylo4(x, check.node.labels="drop", annote=annote)
        res <- phylo4d(tree, tip.data=tip.data, node.data=node.data,
                       all.data=all.data, metadata=metadata, ...)
        res <- addData(res, node.data=nlab.data, pos="before", match.data=FALSE)
    }
    else {
        tree <- phylo4(x, check.node.labels=check.node.labels, annote=annote)
        res <- phylo4d(tree, tip.data=tip.data, node.data=node.data,
                       all.data=all.data, metadata=metadata, ...)
    }

    return(res)
})

### first arg is a phylo4d
##' @rdname phylo4d-methods
##' @aliases phylo4d,phylo4d,phylo4d-method
setMethod("phylo4d", c("phylo4d"), function(x, ...) {
          stop("Your object is already a phylo4d object. If you want to modify",
               " the data attached to it look at the help for tdata()<-,")
      })

### first arg is nexml
##' @rdname phylo4d-methods
##' @aliases nexml,phylo4d-method
setMethod("phylo4d", c("nexml"), function(x) {
    tr <- RNeXML::get_trees_list(x)
    chr <- RNeXML::get_characters(x)
    if (is.null(tr[[1]])) {
        new("phylo4d")
    } else {
        if (length(tr) > 1) {
            warning("Only the first tree has been imported.")
        }
        phylo4d(x=tr[[1]][[1]], chr)
    }
})


### Core function that takes care of the data
.phylo4Data <- function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                        merge.data=TRUE) {

    ## Check validity of phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## Create placeholder data frames for any null data arguments
    if (is.null(tip.data)) tip.data <- formatData(x, NULL, "tip")
    if (is.null(node.data)) node.data <- formatData(x, NULL, "internal")
    if (is.null(all.data)) all.data <- formatData(x, NULL, "all")

    # don't allow all.data columns of same name as tip.data or node.data
    colnamesTipOrNode <- union(names(tip.data), names(node.data))
    if (any(names(all.data) %in% colnamesTipOrNode)) {
        stop("all.data column names must be distinct from ",
             "tip.data and node.data column names")
    }

    ## combine common columns and move into all.data if merging,
    ## otherwise rename them
    colsToMerge <- intersect(names(tip.data), names(node.data))
    if (merge.data && length(colsToMerge)>0) {
        ##TODO could really just index rows directly on 1:nTip and
        ## (nTip+1):(nTip+nNode) in the next two statements for speed,
        ## but this is more robust to changes in node numbering rules
        tip.rows <- tip.data[match(nodeId(x, "tip"),
            row.names(tip.data)), colsToMerge, drop=FALSE]
        node.rows <- node.data[match(nodeId(x, "internal"),
            row.names(tip.data)), colsToMerge, drop=FALSE]
        merge.data <- rbind(tip.rows, node.rows)
        all.data <- data.frame(all.data, merge.data)
    } else {
        names(tip.data)[names(tip.data) %in% colsToMerge] <-
            paste(colsToMerge, "tip", sep=".")
        names(node.data)[names(node.data) %in% colsToMerge] <-
            paste(colsToMerge, "node", sep=".")
    }
    ## now separate tips-only and nodes-only data
    tip.only.data <- tip.data[setdiff(names(tip.data), names(node.data))]
    node.only.data <- node.data[setdiff(names(node.data), names(tip.data))]

    ## combine all data
    complete.data <- data.frame(all.data, tip.only.data, node.only.data)

    ## drop any rows that only contain NAs
    if (ncol(complete.data)==0) {
        return(data.frame())
    } else {
        empty.rows <- as.logical(rowSums(!is.na(complete.data)))
        return(complete.data[empty.rows, , drop=FALSE])
    }

}
