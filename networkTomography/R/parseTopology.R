
#' Build routing matrix from table of link relationships
#'
#' Constructs routing matrix from link relationships. Determines routes using
#' (weighted) shortest-path calculation (mirroring OSPF). Currently handles tied
#' paths arbitrarily; will incorporate fractions for tie resolution in next
#' version. Can optionally include aggregate source and destination flows for
#' each node; this can make a major difference for some topologies. Tomogravity
#' methods typically make use of such information, which most routers collect.
#' Note that resulting routing matrix need not be of full row rank.
#'
#' @param nodes vector (lenght n) of node identifiers
#' @param src vector (length m) of sources, one per link, matched with dest
#' @param dest vector (length n) of destination identifiers, one per link,
#'          matched with src
#' @param weights numeric vector (length m) of weights for each link; used in
#'          shortest-path routing calculations (roughly OSPF)
#' @param agg logical for whether to include aggregate source and destination
#'          flows for each node
#' @param sep character separator between node id's for link and OD names
#' @param aggChar character to indicate aggregate flows; should be distinct from
#'          sep
#' @param verbose integer level of verbosity; 0 is silent, >=1 are increasing
#'          levels of reporting
#' @return List consisting of routing matrix \code{A} (dense) of dimensions m x
#'        n and iGraph object for network \code{topo}
#' @keywords array
#' @export
buildRoutingMatrix <- function(nodes, src, dest, weights=NULL, agg=FALSE,
                               sep='_', aggChar='*', verbose=0) {
    # Format input
    nodes <- as.character(nodes)
    src <- as.character(src)
    dest <- as.character(dest)

    # Get dimensions
    n <- length(nodes)
    m <- length(src)

    # Setup routing matrix
    A <- matrix(integer(1), m+agg*2*n, n*n)

    linkNames <- paste(src, dest, sep=sep)

    if (agg) {
        # Aggregate in flows
        linkNames <- c(linkNames, paste(aggChar, nodes, sep=sep))

        # Aggregate out flows
        linkNames <- c(linkNames, paste(nodes, aggChar, sep=sep))
    }

    rownames(A) <- linkNames
    
    od <- expand.grid(nodes, nodes)
    odNames <- apply(od, 1, paste, sep=sep, collapse=sep)
    colnames(A) <- odNames

    # Build topology (adjacency) matrix
    topo <- matrix(0, n, n)
    colnames(topo) <- rownames(topo) <- nodes

    for (link in 1:m) {
        topo[src[link],dest[link]] <- 1 + topo[src[link],dest[link]]
    }

    # Put topology into igraph object
    topo <- graph.adjacency(topo)
    
    # Iterate over origins
    for (orig in 1:n) {
        # Get list of shortest (weighted) paths
        pathList <- get.shortest.paths(topo, from=nodes[orig], to=nodes,
                                       mode="out", weights=weights)

        # Iterate through destinations
        for (dest in seq(n)[-orig]) {
            # Get links that this OD flow hits
            odName <- paste(nodes[orig], nodes[dest], sep=sep)
            nLinks <- length(pathList[[dest]]) - 1

            links <- sapply(seq(nLinks), function(j)
                            pathList[[dest]][j+c(0,1)] )

            if (verbose > 0)
                cat(sprintf("OD flow: %s\n", odName))

            # Add ones to the correct entries in the corresponding column of A
            for (l in 1:nLinks) {
                linkName <- paste(nodes[links[,l]], collapse=sep)
                
                if (verbose > 1)
                    cat(sprintf("\tLink: %s\tOD: %s\n", linkName, odName),
                        file=stderr())
                A[linkName, odName] <- 1
            }
        }
    }

    # Add in and out flows if requested
    if (agg) {
        # Iterate over nodes
        for (node in nodes) {
            # Get names of in and out flows
            inName <- paste(aggChar, node, sep=sep)
            outName <- paste(node, aggChar, sep=sep)

            # Find in and out OD flows
            inFlows <- grep(paste0(sep, node), colnames(A), fixed=TRUE)
            outFlows <- grep(paste0(node, sep), colnames(A), fixed=TRUE)

            # Add ones to appropriate entries in rows
            A[inName, inFlows] <- 1
            A[outName, outFlows] <- 1
        }
    }

    return(list(A=A, topo=topo))
}
