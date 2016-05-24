## REQUIRED for all trees


##' Validity checking for phylo4 objects
##'
##' Basic checks on the validity of S4 phylogenetic objects
##'
##'
##' @aliases checkPhylo4 checkTree checkPhylo4Data
##' @param object A prospective phylo4 or phylo4d object
##' @return As required by \code{\link[methods]{validObject}}, returns an error
##' string (describing problems) or TRUE if everything is OK.
##' @note
##'
##' These functions are only intended to be called by other phylobase functions.
##'
##' \code{checkPhylo4} is an (inflexible) wrapper for
##' \code{checkTree}.  The rules for \code{phylo4} objects essentially
##' follow those for \code{phylo} objects from the \code{ape} package,
##' which are in turn defined in
##' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
##' These are essentially that: \itemize{ \item if the tree has edge
##' lengths defined, the number of edge lengths must match the number
##' of edges; \item the number of tip labels must match the number of
##' tips; \item in a tree with \code{ntips} tips and \code{nnodes}
##' (total) nodes, nodes 1 to \code{ntips} must be tips \item if the
##' tree is rooted, the root must be node number \code{ntips+1} and
##' the root node must be the first row of the edge matrix \item tip
##' labels, node labels, edge labels, edge lengths must have proper
##' internal names (i.e.  internal names that match the node numbers
##' they document) \item tip and node labels must be unique }
##'
##' You can alter some of the default options by using the function
##' \code{phylobase.options}.
##'
##' For \code{phylo4d} objects, \code{checkTree} also calls
##' \code{checkPhylo4Data} to check the validity of the data associated with the
##' tree. It ensures that (1) the data associated with the tree have the correct
##' dimensions, (2) that the row names for the data are correct.
##' @author Ben Bolker, Steven Kembel, Francois Michonneau
##' @seealso the \code{\link{phylo4}} constructor and
##' \linkS4class{phylo4} class; the \code{\link{phylo4d-methods}} constructor
##' and the \linkS4class{phylo4d} class do checks for the data
##' associated with trees.  See \code{\link{coerce-methods}} for
##' translation functions and \code{\link{phylobase.options} to change
##' some of the default options of the validator.}
##' @include RcppExports.R
##' @include phylo4-class.R
##' @include phylo4-methods.R
##' @export
##' @keywords misc
checkPhylo4 <- function(object) {
    ct <- checkTree(object)

    if (class(object) == "phylo4d")
        ## checkPhyo4Data returns TRUE or fail
        cd <- checkPhylo4Data(object)

    return(ct)
}

checkTree <- function(object) {

    ## case of empty phylo4 object
    if(nrow(object@edge) == 0 && length(object@edge.length) == 0 &&
       length(object@label) == 0 && length(object@edge.label) == 0)
        return(TRUE)

    ## get options
    opt <- phylobase.options()

    ## Storage of error/warning messages
    err <- wrn <- character(0)

    ## Matrix is integer
    if (!is.integer(object@edge)) {
        err <- c(err, "Edge matrix needs to be integer.")
    }

    ## Matrix doesn't have NAs
    if (any(is.na(object@edge))) {
        err <- c(err, "Edge matrix cannot have NAs at this time.",
                 "This could only happen if singletons were allowed",
                 "but this is not supported by phylobase yet.")
    }

    ## Having non-integer or NAs cause cryptic messages, so stop here
    ##  if it's the case
    if (length(err)) return(err)

    ## Named slots
    if (is.null(attributes(object@label)$names)) {
        err <- c(err, "The label slot needs to be a named vector.")
        attributes(object@label) <- list(names=character(0))
    }
    if (is.null(attributes(object@edge.length)$names)) {
        err <- c(err, "The edge.length slot needs to be a named vector.")
        attributes(object@edge.length) <- list(names=character(0))
    }
    if (is.null(attributes(object@edge.label)$names)) {
        err <- c(err, "The edge.label slot needs to be a named vector.")
        attributes(object@edge.label) <- list(names=character(0))
    }

    res <- checkTreeCpp(object, opts=opt)

    if (hasRetic(object)) {
        msg <- "Tree is reticulated."
        if (identical(opt$retic, "fail")) {
            err <- c(err, msg)
        }
        if (identical(opt$retic, "warn")) {
            wrn <- c(wrn, msg)
        }
    }

    if (hasEdgeLength(object) && any(is.na(edgeLength(object)))) {
        naElen <- names(which(is.na(object@edge.length)))
        if (! identical(naElen, edgeId(object, "root")))
            err <- c(err, "Only the root can have NA as edge length. ")
    }

    if (!object@order %in% phylo4_orderings) {
      err <- c(err, paste("unknown order: allowed values are",
               paste(phylo4_orderings,collapse=",")))
    }

    err <- ifelse(nzchar(res[[1]]), c(err, res[[1]]), err)
    wrn <- ifelse(nzchar(res[[2]]), c(wrn, res[[2]]), wrn)

    if (!is.na(wrn)) {
        wrn <- paste(wrn, collapse=", ")
        warning(wrn)
    }
    if (!is.na(err)) {
        err <- paste(err, collapse=", ")
        return(err) #failures are returned as text
    }
    else {
        return(TRUE)
    }

}

## checkTreeOld <- function(object) {

##     ## case of empty phylo4 object
##     if(nrow(object@edge) == 0 && length(object@edge.length) == 0 &&
##        length(object@label) == 0 && length(object@edge.label) == 0)
##         return(TRUE)

##     ## get options
##     opt <- phylobase.options()

##     ## Storage of error/warning messages
##     err <- wrn <- character(0)

##     ## Define variables
##     nedges <- nEdges(object)
##     ntips <- nTips(object)
##     E <- edges(object)
##     tips <- unique(sort(E[,2][!E[,2] %in% E[,1]]))
##     nodes <- unique(sort(c(E)))
##     intnodes <- nodes[!nodes %in% tips]
##     nRoots <- length(which(E[,1] == 0))

##     ## Check edge lengths
##     if (hasEdgeLength(object)) {
##         if (length(object@edge.length) != nedges)
##             err <- c(err, "edge lengths do not match number of edges")
##         ##if(!is.numeric(object@edge.length)) # not needed
##         ##  err <- c(err, "edge lengths are not numeric")
##         ## presumably we shouldn't allow NAs mixed
##         ## with numeric branch lengths except at the root
##         if (sum(is.na(object@edge.length)) > (nRoots + 1))
##             err <- c(err, "NAs in edge lengths")
##         ## Strip root edge branch length (if set to NA)
##         if (any(object@edge.length[!is.na(object@edge.length)] < 0))
##             err <- c(err, "edge lengths must be non-negative")
##         ## Check edge length labels
##         elen.msg <- "Use edgeLength<- to update them."
##         if (is.null(names(object@edge.length))) {
##             err <- c(err, paste("Edge lengths must have names matching edge IDs.",
##                                 elen.msg))
##         }
##         if (!all(names(object@edge.length) %in% edgeId(object, "all"))) {
##             err <- c(err, paste("One or more edge lengths has an unmatched ID name.",
##                                 elen.msg))
##         }
##     }

##     ## Make sure tips and
##     if (!(all(tips==1:ntips) && all(nodes=(ntips+1):(ntips+length(intnodes)))))
##         err <- c(err, "tips and nodes incorrectly numbered")

##     ##careful - nAncest does not work for counting nRoots in unrooted trees
##     nAncest <- tabulate(na.omit(E)[, 2],nbins=max(nodes)) ## bug fix from Jim Regetz
##     nDesc <- tabulate(na.omit(E[,1]))
##     nTips <- sum(nDesc==0)
##     if (!all(nDesc[1:nTips]==0))
##         err <- c(err, "nodes 1 to nTips must all be tips")

##     if (nRoots > 0) {
##       if (sum(E[, 1] == 0) != 1) {
##           err <- c(err, "for a rooted tree, edge matrix must contain (exactly one) explicit root edge with ancestor==0")
##       }
##       root.node <- unname(E[which(E[,1] == 0), 2])
##     }

##     ## Check that nodes are correctly numbered
##     if (!all(nDesc[(nTips+1):(nTips+nNodes(object))]>0))
##         err <- c(err, "nodes (nTips+1) to (nTips+nNodes) must all be internal nodes")

##     ## how do we identify loops???
##     ## EXPERIMENTAL: could be time-consuming for large trees?
##     if (FALSE) {
##       Emat <- matrix(0,nrow=max(E),ncol=max(E))
##       Emat[E] <- 1
##     }
##     if (!object@order %in% phylo4_orderings) {
##       err <- c(err, paste("unknown order: allowed values are",
##                paste(phylo4_orderings,collapse=",")))
##     }

##     ## make sure tip/node labels have internal names that match node IDs
##     lab.msg <- "Use tipLabels<- (and nodeLabels<- if needed) to update them."
##     if (is.null(names(object@label))) {
##         err <- c(err, paste("Tip and node labels must have names matching node IDs.",
##                             lab.msg))

##     } else {
##         if (!all(tips %in% names(na.omit(object@label)))) {
##             err <- c(err, paste("All tips must have associated tip labels.",
##                                 lab.msg))
##         }
##         if (!all(names(object@label) %in% nodeId(object, "all"))) {
##             err <- c(err, paste("One or more tip/node label has an unmatched ID name",
##                                 lab.msg))
##         }
##     }

##     ## make sure edge labels have internal names that match the edges
##     elab.msg <- "Use edgeLabels<- to update them."
##     if(hasEdgeLabels(object)) {
##         if (is.null(names(object@edge.label))) {
##             err <- c(err, paste("Edge labels must have names matching edge IDs.",
##                                 elab.msg))
##         }
##         if (!all(names(object@edge.label) %in% edgeId(object, "all"))) {
##             err <- c(err, paste("One or more edge labels has an unmatched ID name.",
##                                 elab.msg))
##         }
##     }

##     ## make sure that tip and node labels are unique
##     if (hasDuplicatedLabels(object)) {
##         currmsg <- "Labels are not unique"
##         if (opt$allow.duplicated.labels == "fail")
##             err <- c(err, currmsg)
##         if (opt$allow.duplicated.labels == "warn")
##             wrn <- c(wrn, currmsg)
##     }

##     if (any(nDesc>2)) {
##         currmsg <- "tree includes polytomies"
##         if (opt$poly == "fail")
##             err <- c(err, currmsg)
##         if (opt$poly == "warn")
##             wrn <- c(wrn, currmsg)
##       }

##     if (nRoots>1) {
##         currmsg <- "tree has more than one root"
##         if (opt$multiroot == "fail")
##             err <- c(err, currmsg)
##         if (opt$multiroot == "warn")
##             wrn <- c(wrn,currmsg)
##     }
##     if (any(nDesc==1)) {
##         currmsg <- "tree contains singleton nodes"
##         if (opt$singleton == "fail")
##             err <- c(err, currmsg)
##         if (opt$singleton == "warn")
##             wrn <- c(wrn, currmsg)
##     }
##     if (any(nAncest>1)) {
##       currmsg <- paste("tree is reticulated [most functions in phylobase haven't",
##                        "been tested with reticulated trees]")
##       if (opt$retic == "fail")
##           err <- c(err, currmsg)
##       if (opt$retic == "warn")
##           wrn <- c(wrn, currmsg)
##     }
##     if (length(wrn) > 0) {
##         wrn <- paste(wrn, collapse=", ")
##         warning(wrn)
##     }
##     if (length(err) > 0) {
##         err <- paste(err, collapse=", ")
##         return(err) #failures are returned as text
##     }
##     else {
##         return(TRUE)
##     }
## }

checkPhylo4Data <- function(object) {

    ## These are just some basic tests to make sure that the user does not
    ## alter the object in a significant way

    ## Check rownames
    if (nrow(object@data) > 0 &&
        !all(row.names(object@data) %in% nodeId(object, "all")))
        stop("The row names of tree data do not match the node numbers")

    return(TRUE)
}
