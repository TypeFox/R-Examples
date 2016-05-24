#' Convert time string to decimal hour
#'
#' @param x input character vector of times
#' @param fmt input character format for times
#' @return numeric vector of decimal times in hours
#' @keywords character
#' @export
#' @examples
#' strphour("31/08/87 12:53:29")
#'
strphour <- function(x, fmt="(%m/%d/%y %H:%M:%S)") {
    ptime <- strptime(as.character(x), format=fmt)
    ptime$hour + ptime$min/60 + ptime$sec/3600
}

#' Build routing matrix for star network topology
#' 
#' @param n integer number of nodes in the network
#' @return matrix of dimension 2n x n^2 that transforms OD flows to link loads
#' @keywords array
#' @export
#' @examples
#' buildStarMat(3)
#'
buildStarMat <- function(n) {
    # Allocate matrix
    p <- n*n
    J <- n*2
    routemat <- matrix(0, J, p)

    # Setup nonzero entries
    for (i in 1:(nrow(routemat)/2)) {
        routemat[i,seq((i-1)*J/2+1, i*J/2)] <- 1
        routemat[i+nrow(routemat)/2,seq(i,ncol(routemat),J/2)] <- 1
    }

    return(routemat)
}

#' Build routing matrices for linked star topologies; that is, a set of
#' star-topology networks with links between a subset of routers
#'
#' @param nVec integer vector containing number of nodes in each sub-network
#'      (length m)
#' @param Cmat matrix (m x m) containing a one for each linked sub-network;
#'      only upper triangular part is used
#' @return routing matrix of dimension at least 2*sum(nVec) x sum(nVec^2)
#' @seealso \code{\link{buildStarMat}}, which this function depends upon
#' @export
#' @examples
#' nVec <- c(3, 3, 3)
#' Cmat <- diag(3)
#' Cmat[1,2] <- Cmat[2,3] <- 1
#' buildRoutingMat(nVec, Cmat)
buildRoutingMat <- function(nVec, Cmat) {
    # Create overall structure
    N <- sum(nVec)
    m <- length(nVec)

    k <- N*N
    b <- sum(upper.tri(Cmat))
    l <- 2*(N+b)

    A <- matrix(0, l, k)

    # Build star-topology routing matrix for all nodes
    A[1:(2*N),] <- buildStarMat(N)

    # Build between-router part of routing matrix
    n <- 1
    for (i in 1:(m-1)) {
        for (j in (i+1):m) {
            if (Cmat[i,j] == 1) {
                # Create row for i -> j
                srcStart <- c(0,cumsum(nVec))[i] + 1
                srcEnd <- cumsum(nVec)[i]

                dstStart <- c(0,cumsum(nVec))[j] + 1 + N
                dstEnd <- cumsum(nVec)[j] + N

                A[ 2*N + n, ] <- (colSums(A[srcStart:srcEnd,]) * 
                                  colSums(A[dstStart:dstEnd,]) )

                # Create row for j -> i
                srcStart <- c(0,cumsum(nVec))[j] + 1
                srcEnd <- cumsum(nVec)[j]

                dstStart <- c(0,cumsum(nVec))[i] + 1 + N
                dstEnd <- cumsum(nVec)[i] + N

                A[ 2*N + n + 1, ] <- (colSums(A[srcStart:srcEnd,]) * 
                                      colSums(A[dstStart:dstEnd,]) )

                # Augment n
                n <- n + 2
            }
        }
    }

    return(A)
}

#' Make diagonal matrix from vector
#'
#' Build matrix with supplied vector on diagonal; this is much faster than diag
#' due to the use of matrix instead of array
#'
#' @param x numeric vector for diagonal
#' @return matrix of size length(x) x length(x) with x along diagonal
#' @keywords array
#' @seealso \code{\link{diag_ind}}
#' @export
#' @examples
#' diag_mat(seq(5))
diag_mat <- function(x) {
    n <- length(x)
    y <- matrix(0,n,n)
    y[1L + 0L:(n-1L) * (n+1L)] <- x
    return(y)
}

#' Make vector of 1-dimensional diagonal indices for square matrix
#'
#' Compute vector of indices for efficient access to diagonal of a square matrix
#'
#' @param n integer dimension of (square) matrix
#' @return integer vector of length n with indices (unidimensional) of square
#'      matrix
#' @keywords array
#' @seealso \code{\link{diag_mat}}
#' @export
#' @examples
#' ind <- diag_ind(5)
#' diag_mat(seq(5))[ind]
diag_ind <- function(n) {
    1L + 0L:(n-1L)*(n+1L)
}

#' Thinning vector of indices for MCMC
#'
#' Returns a vector of indices with a given spacing for thinning MCMC results
#'
#' @param m integer length of results
#' @param interval thinning interval
#' @return integer vector of indices for thinning
#' @keywords manip ts
#' @export
thin <- function(m, interval=10) {
    seq(1,m,interval)
}

#' Function to aggregate results from matrix to matrix
#'
#' Defaults to mean, SD, limits, and given quantiles. Used to limit memory
#' consumption from MCMC runs.
#'
#' @param mat input numeric matrix to summarize
#' @param q quantiles of mat's columns to provide in summary matrix
#' @return matrix with each row corresponding to a summary measure and each
#'      column corresponding to a column of mat
#' @keywords manip arith
#' @export
#' @examples
#' mat <- matrix(rnorm(5e3), ncol=5)
#' agg(mat)
agg <- function(mat, q=c(0.05, 0.16, 0.5, 0.84, 0.95)) {
    # Convert to matrix if needed
    if (is.vector(mat)) {
        mat <- as.matrix(mat)
    }

    ans <- matrix(0, length(q)+4, ncol(mat))
    nc <- ncol(mat)

    # Dimnames
    dimnames(ans) <- list()
    dimnames(ans)[[1]] <- vector("character", nrow(ans))
    dimnames(ans)[[1]][1:4] <- c("mean", "sd", "min", "max")

    # Default aggregates
    ans[1,] <- colMeans(mat)
    ans[2,] <- apply(mat, 2, sd)
    ans[3:4,] <- sapply( 1:nc, function(j) range(mat[,j]) )

    # Quantiles
    if (length(q) > 0) {
        ans[5:nrow(ans),] <- sapply(1:nc, function(j)
                                    quantile(mat[,j], q) )
        dimnames(ans)[[1]][5:nrow(ans)] <- sprintf("q%02d", q*100)
    }

    return(ans)
}

#' Check for deterministically-known OD flows at single time
#'
#' Uses xranges from limSolve to find deterministically-known OD flows
#'
#' @param y numeric vector of link loads, dimension m
#' @param A routing matrix of dimension m x k
#' @return logical vector of length k; TRUE for unknown OD flows, FALSE for
#'      known
#' @keywords algebra
#' @export
#' @examples
#' data(bell.labs)
#' getActive(bell.labs$Y[1,], bell.labs$A)
getActive <- function(y, A) {
    odRanges <- xranges(E=A, F=y, ispos=TRUE)
    activeOD <- (odRanges[,2]-odRanges[,1]>0)
    return( activeOD )
}

#' Find indices of source and destination for each point-to-point flow
#' 
#' This works only for routing matrices that include all aggregate source and
#' destination flows. It is often easier to build these indices manually via
#' string processing or during the construction of the routing matrix.
#' 
#' @param A routing matrix of dimension m x k. This should be the reduced-rank 
#'   version including all aggregate source and destination flows.
#' @return list consisting of two component, src and dst, which are integer 
#'   vectors of length k containing the index (in y = A x) of the source and 
#'   destination flows that each point-to-point flow is part of.
#' @keywords algebra
#' @export
#' @examples
#' data(cmu)
#' src.dst.ind <- getSrcDstIndices(cmu$A.full)
getSrcDstIndices <- function(A) {
    indList <- list(src=integer(ncol(A)), dst=integer(ncol(A)))

    for (i in 1:ncol(A)) {
        # Get link loads that include xi
        ind <- which(A[,i] > 0)
        
        # Get number of flows in each of these loads; assuming these are source
        # and destination
        l <- sapply(ind, function(j) sum(A[j,] > 0))

        srcDestInd <- ind[which(rank(l, ties.method="min") <= 2)]
        indList$src[i] <- srcDestInd[1]
        indList$dst[i] <- srcDestInd[2]
    }
    return(indList)
}

#' Compute pivoted decomposition of routing matrix A into full-rank and
#' remainder, as in Cao et al. 2000, via the QR decomposition.
#' 
#' @param A routing matrix of dimension m x k
#' @return list containing two matrices: A1 (m x m), a full-rank subset of the
#'      columns of A, and A2 (m x k - m), the remaining columns
#' @keywords algebra
#' @export
decomposeA <- function(A) {
    Aqr <- qr(A)
    A1  <- A[,Aqr$pivot[seq(Aqr$rank)]]
    A2  <- A[,Aqr$pivot[seq(Aqr$rank+1,ncol(A))]]
    return(list(A1=A1, A2=A2))
}

#' Compute total traffic from a particular time.
#' 
#' @param yt length-m numeric vectors of observed aggregate flows at a
#'      particular time
#' @param A1 m x m matrix containing the full-rank portion of the network's
#'      routing matrix, as supplied by \code{\link{decomposeA}}
#' @keywords algebra
#' @export
#' @examples
#' data(bell.labs)
#' A.decomp <- decomposeA(bell.labs$A)
#' total.traffic <- calcN(yt=bell.labs$Y[1,], A1=A.decomp$A1)
#' total.traffic == sum(bell.labs$X[1,])
calcN <- function(yt, A1) {
    sum(qr.solve(A1, yt))
}
