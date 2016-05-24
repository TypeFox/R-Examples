##
## INPUT:
## boot: number of bootstrap iterations
## report: whether to print status bar
## estimate: original fit coefficients
## y: dependent variable
## a: acceptance vector (for ultimatum only)
## regr: list of regressor matrices
## fn: log-likelihood function
## gr: gradient function (if any)
## fixed: logical vector indicating which parameters are held fixed
## method: optimization routine to use
##
## RETURN:
## matrix of bootstrap results, each row an iteration
## 
gameBoot <- function(boot, report = TRUE, estimate, y, a = NULL, regr, fn, gr,
                     fixed, method, ...)
{
    bootMatrix <- matrix(NA, nrow = boot, ncol = length(estimate))
    failedBoot <- logical(boot)
    if (report) {
        cat("\nRunning bootstrap iterations...\n")
        pb <- txtProgressBar(min = 1, max = boot)
    }
    for (i in seq_len(boot)) {
        bootSamp <- sample(seq_len(length(y)), replace = TRUE)
        newy <- y[bootSamp]
        newa <- a[bootSamp]  ## for the ultimatum model
        newregr <- lapply(regr, function(x) x[bootSamp, , drop = FALSE])
        bootResults <- maxLik(logLik = fn, grad = gr, start = estimate, fixed =
                              fixed, method = method, y = newy, acc = newa, regr
                              = newregr, ...)
        cc <- convergenceCriterion(method)
        if (!(bootResults$code %in% cc)) {
            warning("bootstrap iteration ", i,
                    " failed to converge and will be removed")
            failedBoot[i] <- TRUE
        }
        bootMatrix[i, ] <- bootResults$estimate

        if (report)
            setTxtProgressBar(pb, i)
    }
    if (report)
        cat("\n")
    bootMatrix <- bootMatrix[!failedBoot, , drop = FALSE]
    colnames(bootMatrix) <- names(estimate)
    return(bootMatrix)
}
