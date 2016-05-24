
##' Tree traversal and utility functions
##'
##' Functions for describing relationships among phylogenetic nodes (i.e.
##' internal nodes or tips).
##'
##' \code{ancestors} and \code{descendants} can take \code{node} vectors of
##' arbitrary length, returning a list of output vectors if the number of valid
##' input nodes is greater than one. List element names are taken directly from
##' the input node vector.
##'
##' If any supplied nodes are not found in the tree, the behavior currently
##' varies across functions.
##' \itemize{
##' \item Invalid nodes are automatically omitted by \code{ancestors}
##' and \code{descendants}, with a warning.
##'
##' \item \code{ancestor}
##' will return \code{NA} for any invalid nodes, with a warning.
##'
##' \item Both \code{children} and \code{siblings} will return an empty
##' vector, again with a warning.
##' }
##' @param phy a \linkS4class{phylo4} object (or one inheriting from
##' \linkS4class{phylo4}, e.g. a \linkS4class{phylo4d} object)
##' @param node either an integer corresponding to a node ID number, or a
##' character corresponding to a node label; for \code{ancestors} and
##' \code{descendants}, this may be a vector of multiple node numbers or names
##' @param type (\code{ancestors}) specify whether to return just direct
##' ancestor ("parent"), all ancestor nodes ("all"), or all ancestor nodes
##' including self ("ALL"); (\code{descendants}) specify whether to return just
##' direct descendants ("children"), all extant descendants ("tips"), or all
##' descendant nodes ("all") or all descendant nodes including self ("ALL").
##' @param include.self whether to include self in list of siblings
##' @param \dots a list of node numbers or names, or a vector of node numbers or
##' names
##' @return \describe{
##' \item{\code{ancestors}}{ return a named vector (or a list
##' of such vectors in the case of multiple input nodes) of the
##' ancestors and descendants of a node}
##'
##' \item{\code{descendants}}{ return a named vector (or a list of
##' such vectors in the case of multiple input nodes) of the ancestors
##' and descendants of a node}
##'
##' \item{\code{ancestor}}{ \code{ancestor} is analogous to
##' \code{ancestors(\dots{}, type="parent")} (i.e. direct ancestor
##' only), but returns a single concatenated vector in the case of
##' multiple input nodes}
##'
##' \item{\code{children}}{is analogous to \code{descendants(\dots{},
##' type="children")} (i.e. direct descendants only), but is not
##' currently intended to be used with multiple input nodes }
##'
##' \item{\code{siblings}}{ returns sibling nodes (children of the same
##' parent)}
##' }
##' @seealso \code{\link[ape]{mrca}}, in the ape package, gives a list of all
##' subtrees
##' @export
##' @rdname ancestors
##' @include phylo4-class.R
##' @include phylo4-methods.R
##' @include getNode-methods.R
##' @examples
##'
##'   data(geospiza)
##'   nodeLabels(geospiza) <- LETTERS[1:nNodes(geospiza)]
##'   plot(as(geospiza, "phylo4"), show.node.label=TRUE)
##'   ancestor(geospiza, "E")
##'   children(geospiza, "C")
##'   descendants(geospiza, "D", type="tips")
##'   descendants(geospiza, "D", type="all")
##'   ancestors(geospiza, "D")
##'   MRCA(geospiza, "conirostris", "difficilis", "fuliginosa")
##'   MRCA(geospiza, "olivacea", "conirostris")
##'
##'   ## shortest path between 2 nodes
##'   shortestPath(geospiza, "fortis", "fuliginosa")
##'   shortestPath(geospiza, "F", "L")
##'
##'   ## branch length from a tip to the root
##'   sumEdgeLength(geospiza, ancestors(geospiza, "fortis", type="ALL"))
ancestor <- function(phy,node) {
    node2 <- getNode(phy,node)
    ## r <- which(edges(phy)[,2]==node)
    r <- match(node2,edges(phy)[,2])
    return(getNode(phy,edges(phy)[r,1],missing="OK"))
}


##' @rdname ancestors
##' @aliases children
##' @export
children <- function(phy,node) {
    node2 <- getNode(phy,node)
    r <- which(edges(phy)[,1]==node2)
    getNode(phy,edges(phy)[r,2])
}

##' @rdname ancestors
##' @aliases descendants
##' @export
descendants <- function (phy, node, type=c("tips","children","all", "ALL")) {
    type <- match.arg(type)

    ## look up nodes, warning about and excluding invalid nodes
    oNode <- node
    node <- getNode(phy, node, missing="warn")
    isValid <- !is.na(node)
    node <- as.integer(node[isValid])

    if (type == "children") {
        res <- lapply(node, function(x) children(phy, x))
        names(res) <- node
    } else {
        ## edge matrix must be in preorder for the C function!
        if (phy@order=="preorder") {
            edge <- phy@edge
        } else {
            edge <- reorder(phy, order="preorder")@edge
        }
        ## extract edge columns
        ancestor <- as.integer(edge[, 1])
        descendant <- as.integer(edge[, 2])

        ## return indicator matrix of ALL descendants (including self)
        isDes <- .Call("descendants", node, ancestor, descendant)
        storage.mode(isDes) <- "logical"

        if (type == "all") {
            i <- match(intersect(node, nodeId(phy, "internal")), descendant)
            isDes[i, seq_along(node)] <- FALSE
        }

        ## if only tips desired, drop internal nodes
        if (type=="tips") {
            isDes[descendant %in% nodeId(phy, "internal"),] <- FALSE
        }
        res <- lapply(seq_along(node), function(n) {
                          getNode(phy, descendant[isDes[,n]])
                      })
        names(res) <- node
    }

    ## if just a single node, return as a single vector
    if (length(res)==1) res <- res[[1]]
    res

    ## Original pure R implementation of the above
    ## (note that it does not require preorder ordering)
    ##n <- nTips(phy)
    ##if (node <= n) {
    ##    return(node)
    ##}
    ##l <- numeric()
    ##d <- children(phy, node)
    ##for (j in d) {
    ##    if (j <= n)
    ##      l <- c(l,j)
    ##    else if (type=="all") l <- c(l,j,
    ##               descendants(phy,j,type="all"))
    ##    else l <- c(l, descendants(phy,j,type=type))
    ##}
}

##' @rdname ancestors
##' @aliases siblings
##' @export
siblings <- function(phy, node, include.self=FALSE) {
    v <- children(phy,ancestor(phy,node))
    if (!include.self) v <- v[v!=getNode(phy,node)]
    v
}

##' @rdname ancestors
##' @aliases siblings
##' @export
ancestors <- function (phy, node, type=c("all","parent","ALL")) {

    type <- match.arg(type)

    ## look up nodes, warning about and excluding invalid nodes
    oNode <- node
    node <- getNode(phy, node, missing="warn")
    isValid <- !is.na(node)
    node <- as.integer(node[isValid])

    if (length(node) == 0) {
      return(NA)
    }

    if (type == "parent") {
        res <- lapply(node, function(x) ancestor(phy, x))
    } else {
        ## edge matrix must be in postorder for the C function!
        if (phy@order=="postorder") {
            edge <- phy@edge
        } else {
            edge <- reorder(phy, order="postorder")@edge
        }
        ## extract edge columns
        ancestor <- as.integer(edge[, 1])
        descendant <- as.integer(edge[, 2])

        ## return indicator matrix of ALL ancestors (including self)
        isAnc <- .Call("ancestors", node, ancestor, descendant)
        storage.mode(isAnc) <- "logical"

        ## drop self if needed
        if (type=="all") {
            isAnc[cbind(match(node, descendant), seq_along(node))] <- FALSE
        }
        res <- lapply(seq_along(node), function(n) getNode(phy,
            descendant[isAnc[,n]]))
    }
    names(res) <- as.character(oNode[isValid])

    ## if just a single node, return as a single vector
    if (length(res)==1) res <- res[[1]]
    res

    ## Original pure R implementation of the above
    ## (note that it does not require preorder ordering)
    ##if (node == rootNode(phy))
    ##    return(NULL)
    ##repeat {
    ##    anc <- ancestor(phy, node)
    ##    res <- c(res, anc)
    ##    node <- anc
    ##    if (anc == n + 1)
    ##        break
    ##}
}
