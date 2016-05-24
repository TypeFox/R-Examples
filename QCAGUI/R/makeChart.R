`makeChart` <-
function(primes = "", configs = "", snames = "") {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    prmat <- is.matrix(primes)
    comat <- is.matrix(configs)
    
    if (prmat + comat  == 2) {
        
        if (!(is.numeric(primes) & is.numeric(configs))) {
            cat("\n")
            stop(simpleError("Matrices have to be numeric.\n\n"))
        }
        
        if (any(primes < 0) | any(configs < 0)) {
            cat("\n")
            stop(simpleError("Matrix values have to be non-negative.\n\n"))
        }
        
        if (any(apply(primes, 1, sum) == 0) | any(apply(configs, 1, sum) == 0)) {
            cat("\n")
            stop(simpleError("Matrices have to be specified at implicants level.\n\n"))
        }
        
        primes2 <- matrix(logical(length(primes)), dim(primes))
        primes2[primes > 0] <- TRUE
        
        mtrx <- sapply(seq(nrow(primes)), function(x) {
            apply(configs, 1, function(y) {
                all(primes[x, primes2[x, ]] == y[primes2[x, ]])
            })
        })
            
        if (nrow(configs) == 1) {
            mtrx <- matrix(mtrx)
        }
        else {
            mtrx <- t(mtrx)
        }
        
        collapse = ifelse(all(nchar(colnames(primes)) == 1) & all(nchar(colnames(configs)) == 1), "", "*")
        
        rownames(mtrx) <- QCA::writePrimeimp(primes, collapse = collapse,  uplow = all(primes < 3) | all(configs < 3))
        colnames(mtrx) <- QCA::writePrimeimp(configs, collapse = collapse, uplow = all(primes < 3) | all(configs < 3))
        
        return(mtrx)
        
    }
    else if (prmat + comat  == 0) {
        tconfigs <- translate(configs, snames)
        if (identical(snames, "")) {
            snames <- colnames(tconfigs)
        }
        tprimes <- translate(primes, snames)
        
        mtrx <- matrix(FALSE, nrow=nrow(tprimes), ncol=nrow(tconfigs))
        
        for (i in seq(nrow(mtrx))) {
            for (j in seq(ncol(mtrx))) {
                tp <- tprimes[i, ]
                tc <- tconfigs[j, ]
                
                if (is.element("mv", names(attributes(tprimes)))) {
                    tpmv <- attr(tprimes, "mv")[i, ]
                    tcmv <- attr(tconfigs, "mv")[j, ]
                    mtrx[i, j] <- all(tp[tp >= 0] == tc[tp >= 0]) & all(tpmv[tp >= 0] == tcmv[tp >= 0])
                }
                else {
                    mtrx[i, j] <- all(tp[tp >= 0] == tc[tp >= 0])
                }
            }
        }
        
        colnames(mtrx) <- rownames(tconfigs)
        rownames(mtrx) <- rownames(tprimes)
        return(mtrx)
    }
    else {
        cat("\n")
        stop(simpleError("Both arguments have to be matrices.\n\n"))
    }
}

