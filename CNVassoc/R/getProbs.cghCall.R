getProbs.cghCall <-
function (x)
{
    if (!inherits(x, "cghCall"))
        stop("object must be of class 'cghCall'")
    requireNamespace("Biobase", quietly = TRUE)
    requireNamespace("CGHbase", quietly = TRUE)
    Clone <- Biobase::featureNames(x)
    Chromo <- CGHbase::chromosomes(x)
    BPstart <- CGHbase::bpstart(x)
    BPend <- CGHbase::bpend(x)
    Calls <- CGHbase::calls(x)
    Probloss <- CGHbase::probloss(x)
    Probnorm <- CGHbase::probnorm(x)
    Probgain <- CGHbase::probgain(x)
    colnam <- c(colnames(Probloss), colnames(Probnorm), colnames(Probgain))
    if (is.null(CGHbase::probamp(x))) {
        allprob <- c()
        ncl <- ncol(Probnorm)
        for (i in 1:ncl) {
            Probsall <- cbind(Probloss[, i], Probnorm[, i], Probgain[,
                i])
            colnames(Probsall) <- c(colnam[i], colnam[ncl + i],
                colnam[2 * ncl + i])
            allprob <- cbind(allprob, Probsall)
        }
    }
    else {
        Probamp <- CGHbase::probamp(x)
        colnam <- c(colnam, colnames(Probamp))
        allprob <- c()
        ncl <- ncol(Probnorm)
        for (i in 1:ncl) {
            Probsall <- cbind(Probloss[, i], Probnorm[, i], Probgain[,
                i], Probamp[, i])
            colnames(Probsall) <- c(colnam[i], colnam[ncl + i],
                colnam[2 * ncl + i], colnam[3 * ncl + i])
            allprob <- cbind(allprob, Probsall)
        }
    }
    allprob2 <- round(allprob, 3)
    probs <- data.frame(Clone, Chromo, BPstart, BPend, allprob2)
    probs
}

