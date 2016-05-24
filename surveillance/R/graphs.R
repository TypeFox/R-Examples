################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Functions concerning graphs: neighbourhood order, adjacency matrix
### These are wrappers around functionality from package "spdep" by Roger Bivand
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 666 $
### $Date: 2013-11-08 15:45:36 +0100 (Fre, 08 Nov 2013) $
################################################################################


### Determine the matrix of neighbourhood orders
### given the binary matrix of first-order neighbours.
### Working horse: spdep::nblag()

nbOrder <- function (neighbourhood, maxlag = 1)
{
    if (!requireNamespace("spdep"))
        stop("package ", dQuote("spdep"),
             " is required to determine neighbourhood orders")

    stopifnot(isScalar(maxlag), maxlag > 0)
    checkNeighbourhood(neighbourhood)
    neighbourhood <- neighbourhood == 1           # convert to binary matrix
    nregions <- nrow(neighbourhood)
    maxlag <- as.integer(min(maxlag, nregions-1)) # upper bound of nb order
    
    if (maxlag == 1L) {
        storage.mode(neighbourhood) <- "integer"
        return(neighbourhood)
    }

    ## manually convert to spdep's "nb" class
    ## region.idxs <- seq_len(nregions)
    ## nb <- lapply(region.idxs, function(i) {
    ##     nbs <- which(neighbourhood[i,])
    ##     if (length(nbs) > 0L) nbs else 0L
    ## })
    ## class(nb) <- "nb"

    ## convert first-order neighbourhood to spdep's "nb" class
    nb <- spdep::mat2listw(neighbourhood)$neighbours
    attr(nb, "region.id") <- NULL

    ## compute higher order neighbours using spdep::nblag()
    nb.lags <- spdep::nblag(nb, maxlag=maxlag)

    ## Side note: fast method to determine neighbours _up to_ specific order:
    ## crossprod(neighbourhood) > 0  # up to second order neighbours (+set diag to 0)
    ## (neighbourhood %*% neighbourhood %*% neighbourhood) > 0  # up to order 3
    ## and so on...

    ## convert to a single matrix
    nbmat <- neighbourhood   # logical first-order matrix
    storage.mode(nbmat) <- "numeric"
    for (lag in 2:maxlag) {
        if (any(spdep::card(nb.lags[[lag]]) > 0L)) { # any neighbours of this order
            nbmat.lag <- spdep::nb2mat(nb.lags[[lag]], style="B",
                                       zero.policy=TRUE)
            nbmat <- nbmat + lag * nbmat.lag
        }
    }
    attr(nbmat, "call") <- NULL
    storage.mode(nbmat) <- "integer"

    ## message about maximum neighbour order by region
    maxlagbyrow <- apply(nbmat, 1, max)
    message("Note: range of maximum neighbour order by region is ",
            paste(range(maxlagbyrow), collapse="-"))

    ## Done
    nbmat
}


### Derive adjacency structure from a SpatialPolygons object
### Working horse: spdep::poly2nb

poly2adjmat <- function (SpP, ..., zero.policy = TRUE)
{
    if (!requireNamespace("spdep"))
        stop("package ", dQuote("spdep"),
             " is required to derive adjacencies from SpatialPolygons")
    nb <- spdep::poly2nb(SpP, ...)
    adjmat <- spdep::nb2mat(nb, style="B", zero.policy=zero.policy)
    attr(adjmat, "call") <- NULL
    colnames(adjmat) <- rownames(adjmat)
    adjmat
}
