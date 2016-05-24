### idcoef.R ---
### Author: Na Li <nali@umn.edu>
### Created: 2006/07/10 15:38:27
### Time-stamp: "Mon Jul 10 12:27:28 CDT 2006 (nali@bass.biostat.umn.edu)"
### Version: $Id$

identity.coefs <- function (samples, pedigree)
{
    pedigree <- as.matrix (pedigree)
    if (ncol (pedigree) < 3) {
        stop ("pedigree has to have at least 3 columns")
    }

    if (is.vector (samples)) {
        u <- unique (samples)
        samples <- t (apply (expand.grid (u, u), 1, sort))
        samples <- samples[!duplicated (samples),]
    } else if (!is.matrix (samples) | ncol (samples) != 2) {
        stop ("`samples' has to be a vector or a matrix of two columns")
        u <- unique (as.vector (samples))
    }
    if (!all (u %in% pedigree[,1])) {
        stop ("Some sample IDs not present in the pedigree! ",
              paste (u[!u %in% pedigree[,1]])) 
    }
    samples <- as.matrix (samples)
    storage.mode (samples) <- "integer"
    coefs <- .Call ("compute_idcoefs",
                    as.integer (pedigree[,1]),
                    as.integer (pedigree[,2]),
                    as.integer (pedigree[,3]),
                    t (samples), PACKAGE = "identity")
    cbind (samples, t (coefs))
}
