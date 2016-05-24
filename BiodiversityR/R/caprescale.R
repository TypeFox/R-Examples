`caprescale` <-
function (x, verbose = FALSE) 
{
    if (is.null(x$CCA)) {
        nr <- nrow(x$CA$u)
        nunpeiv <- x$CA$rank
        neiv <- npeiv <- length(x$CA$eig)
        sumpeiv <- sum(x$CA$eig[1:nunpeiv])
        sumeiv <- sum(x$CA$eig)       
    }
    else {
        nr <- nrow(x$CCA$u)
        nceiv <- length(x$CCA$eig)
        nunpeiv <- x$CA$rank
        nuneiv <- length(x$CA$length)
        npeiv <- nceiv + nunpeiv
        neiv <- nceiv + nuneiv
        sumpeiv <- sum(x$CA$eig[1:nunpeiv]) + sum(x$CCA$eig)
        sumeiv <- sum(x$CA$eig) + sum(x$CCA$eig)
    }
    const <- attributes(scores(x))$const
    if (is.null(x$CCA) == F) {
        x$CCA$v <- x$CCA$v * const
        x$CCA$wa <- x$CCA$wa * const
        x$CCA$u <- x$CCA$u * const
        x$CCA$biplot <- x$CCA$biplot * const
        x$CCA$centroids <- x$CCA$centroids * const
    }
    x$CA$v <- x$CA$v * const
    x$CA$u <- x$CA$u * const
    if (verbose == T) {
        distmat <- as.matrix(vegdist(summary(x, axes = npeiv, 
            scaling = 1)$sites, method = "euc"))
        ssdist <- sum((distmat)^2)/(2 * nrow(distmat))
        if (substr(x$inertia, 1, 4) == "mean") {
            sstot <- sumpeiv * (nr - 1)
            sumeiv <- sumeiv * (nr - 1)
        }else {
            sstot <- sumpeiv
            ssdist <- ssdist / (nr-1)         
        }
        cat("SSTot obtained from sum of all eigenvalues:", sumeiv, "\n")
        cat("SSTot obtained from sum of all positive eigenvalues:", sstot, "\n")
        cat("SSTot reflected by distances among site scores on all axes:", 
            ssdist, "\n")
        if (is.null(x$CCA) == F) {
            distmat <- as.matrix(vegdist(summary(x, axes = nceiv, 
                scaling = 1)$constraints, method = "euc"))
            ssdistcf <- sum((distmat)^2)/(2 * nrow(distmat))
            distmat <- as.matrix(vegdist(summary(x, axes = nceiv, 
                scaling = 1)$sites, method = "euc"))
            ssdistcs <- sum((distmat)^2)/(2 * nrow(distmat))
            distmat <- as.matrix(vegdist(summary(x, axes = npeiv, 
                scaling = 1)$sites[, ((nceiv + 1):npeiv)], method = "euc"))
            ssdistus <- sum((distmat)^2)/(2 * nrow(distmat))
            if (substr(x$inertia, 1, 4) == "mean") {
                sstotc <- sum(x$CCA$eig) * (nr - 1)
                sstotu <- sum(x$CA$eig[1:nunpeiv]) * (nr - 1)               
            }
            else {
                sstotc <- sum(x$CCA$eig)
                ssdistcs <- ssdistcs / (nr-1)
                ssdistcf <- ssdistcf / (nr-1)
                sstotu <- sum(x$CA$eig[1:nunpeiv])
                ssdistus <- ssdistus / (nr-1)
            }
            cat("SSExpl obtained from eigenvalues of constrained axes:", 
                sstotc, "\n")
            cat("SSExpl reflected by distances among site scores on constrained axes (scaling 1):", 
                ssdistcs, "\n")
            cat("SSExpl reflected by distances among fitted site scores on constrained axes (scaling 1):", 
                ssdistcf, "\n")
            cat("SSRes obtained from eigenvalues of positive unconstrained axes:", 
                sstotu, "\n")
            cat("SSRes reflected by distances among site scores on positive unconstrained axes (scaling 1):", 
                ssdistus, "\n")
        }
    }
    return(x)
}


