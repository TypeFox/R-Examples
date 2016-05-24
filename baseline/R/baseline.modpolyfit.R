### modpolyfit_LM2003.R:  Implementation of Modified polyfit method; Lieber
### & Mahadevan-Jansen (2003)
### $Id: baseline.modpolyfit.R 90 2007-11-16 15:04:55Z kristl $
### By Bjørn-Helge Mevik and Kristian Hovde Liland

baseline.modpolyfit <- function(spectra, t, degree = 4, tol = 1e-3, rep = 100) {
    ## Removing dimnames can save time and memory with long spectra:
    dimnames(spectra) <- NULL

    np <- dim(spectra)
    baseline  <- matrix(0, np[1], np[2])
    if (missing(t) || (t == FALSE)) t <- 1:np[2]

    ## The scaling of 1 here is so crossprod(polx) == diag(degree + 1),
    ## which makes the linear regressions below very simple:
    polx <- cbind(1/sqrt(np[2]), stats::poly(t, degree = degree))

    for(i in 1:np[1]){
        ywork <- yold <- yorig <- spectra[i,]

        nrep <- 0
        repeat {
            nrep <- nrep + 1
            ## ypred is fitted(lm.fit(x, y)), only faster:
            ypred <- polx %*% crossprod(polx, yold)
            ywork <- pmin(yorig, ypred)
            crit <- sum(abs((ywork - yold) / yold), na.rm = TRUE)
            if (crit < tol || nrep > rep)
                break
            yold <- ywork
        }
        baseline[i,]  <- ypred
    }
    list(baseline = baseline, corrected = spectra - baseline)
}
