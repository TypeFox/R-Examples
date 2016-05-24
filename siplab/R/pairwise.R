pairwise <-
function(plants, maxN=NULL, maxR=NULL, select=NULL, selpar=NULL, kernel, kerpar=NULL) {
# Pairwise competition indices
    # Competition index computation
    cindex <- function(Y, current, dists, dranks,
                       select, selpar, kernel, kerpar) {
        # Narrow down to actual competitors
        if(!is.null(select) && npoints(Y) > 0) {
            comp <- select(imarks=marks(current), jmarks=marks(Y),
                           dists=dists, dranks=dranks, par=selpar)
            Y <- Y[comp]
            dists <- dists[comp]
            dranks <- dranks[comp]
        }
        # Competition index
        if (npoints(Y) > 0)
            return(sum(kernel(imarks=marks(current), jmarks=marks(Y),
                          dists=dists, dranks=dranks, par=kerpar)))
        else return(0)
    } # end of function
    # Get indices
    marks(plants) <- data.frame(marks(plants), cindex=0)
    marks(plants)$cindex <- applynbd(plants, cindex, N=maxN, R=maxR, exclude=TRUE,
            select=select, selpar=selpar, kernel=kernel, kerpar=kerpar)
    return(plants)
}
