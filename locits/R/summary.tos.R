summary.tos <-
function (object, size = 0.05, quiet = FALSE, mctype = "FDR", ...) 
{
    hwtosop <- object
    ntests <- length(hwtosop$spvals)
    bonreject <- sum(hwtosop$spvals < size/ntests)
    rejpval <- hwtosop$rejpval
    if (quiet == FALSE) {
        cat("There are ", ntests, " hypothesis tests altogether\n")
        cat("There were ", hwtosop$nreject, " FDR rejects\n")
        if (hwtosop$nreject == 0) 
            cat("No p-values were smaller than the FDR val of: ", 
                rejpval, "\n")
        else cat("The rejection p-value was ", rejpval, "\n")
        cat("Using Bonferroni rejection p-value is ", size/ntests, 
            "\n")
        cat("And there would be ", bonreject, " rejections.\n")
    }
    ll <- length(hwtosop$AllTS)
    v <- NULL
    count <- 1
    if (hwtosop$nreject != 0) {
        if (mctype == "FDR") {
            local.rejpval <- rejpval
            if (quiet == FALSE) 
                cat("Listing FDR rejects... (thanks Y&Y!)\n")
        }
        else if (mctype == "BON") {
            local.rejpval <- size/ntests
            if (quiet == FALSE) 
                cat("Listing Bonferroni rejects... (thanks B!)\n")
        }
        else stop("Unknown mctype value. Can only be FDR or BON")
        for (i in 1:ll) {
            if (!is.null(hwtosop$AllPVal[[i]])) {
                pvals <- hwtosop$AllPVal[[i]]
                for (k in 0:(pvals$nlevels - 1)) {
                  the.pvals <- accessD(pvals, lev = k)
                  ix <- which(the.pvals <= local.rejpval)
                  if (length(ix) > 0) {
                    if (quiet == FALSE) {
                      cat("P:", i - 1, "HWTlev: ", k, " indices on next line...")
                      print(ix)
                    }
                    v[[count]] <- c(i - 1, k, ix)
                    count <- count + 1
                  }
                }
            }
        }
    }
    if (mctype == "FDR") 
        vret <- list(rejlist = v, nreject = hwtosop$nreject, 
            mctype = mctype)
    else if (mctype == "BON") 
        vret <- list(rejlist = v, nreject = bonreject, mctype = mctype)
    invisible(vret)
}
