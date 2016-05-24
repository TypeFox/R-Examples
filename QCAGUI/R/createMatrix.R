`createMatrix` <-
function(noflevels) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    conds <- length(noflevels)
    pwr <- unique(noflevels)
    
    if (any(pwr > 2)) {
        logical <- FALSE
    }
    
    if (length(pwr) == 1) {
        create <- function(idx) {
            rep.int(c(sapply(seq_len(pwr) - 1, function(x) rep.int(x, pwr^(idx - 1)))),
                    pwr^conds/pwr^idx)
        }
        retmat <- sapply(rev(seq_len(conds)), create)
    }
    else {
        mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
        orep  <- cumprod(rev(c(rev(noflevels)[-1], 1)))
        retmat <- sapply(seq_len(conds), function(x) {
           rep.int(rep.int(seq_len(noflevels[x]) - 1, rep.int(mbase[x], noflevels[x])), orep[x])
        })
    }
    
    if (is.vector(retmat)) {
        retmat <- matrix(retmat, nrow=1)
    }
    
    return(retmat)
}

