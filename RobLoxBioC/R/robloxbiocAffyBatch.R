###############################################################################
## Use robloxbioc to preprocess Affymetrix-data - comparable to MAS 5.0
###############################################################################
setMethod("robloxbioc", signature(x = "AffyBatch"),
    function(x, bg.correct = TRUE, pmcorrect = TRUE, normalize = FALSE,
            add.constant = 32, verbose = TRUE, 
            eps = NULL, eps.lower = 0, eps.upper = 0.05, steps = 3L, 
            fsCor = TRUE, mad0 = 1e-4, contrast.tau = 0.03, scale.tau = 10, 
            delta = 2^(-20), sc = 500) {
        if(bg.correct){
            if(verbose) cat("Background correcting ...")
            x <- affy::bg.correct.mas(x, griddim = 16)
            if(verbose) cat(" done.\n")
        }
        n <- length(x)
        ids <- featureNames(x)
        m <- length(ids)

        CDFINFO <- getCdfInfo(x)
        INDEX <- sapply(ids, get, envir = CDFINFO)
        NROW <- unlist(lapply(INDEX, nrow))
        nr <- as.integer(names(table(NROW)))

        intensData <- intensity(x)
        rob.est <- matrix(NA, ncol = 2, nrow = m*n)
        if(pmcorrect){
            if(verbose) cat("PM/MM correcting ...")
            diff.log2 <- function(INDEX, x){
                l.pm <- INDEX[,1]
                if(ncol(INDEX) == 2)
                    l.mm <- INDEX[,2]
                else
                    l.mm <- integer()
                log2(x[l.pm, , drop = FALSE]/x[l.mm, , drop = FALSE])
            }
            res <- lapply(INDEX, diff.log2, x = intensData)
            rob.est1 <- matrix(NA, ncol = 2, nrow = m*n)
            for(k in nr){
                ind <- which(NROW == k)
                temp <- matrix(do.call(rbind, res[ind]), nrow = k)
                ind1 <-  as.vector(sapply(seq_len(n)-1, function(x, ind, m){ ind + x*m }, ind = ind, m = m))
                rob.est1[ind1, 1:2] <- robloxbioc(t(temp), eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                                                  steps = steps, fsCor = fsCor, mad0 = mad0)
            }
            sb <- matrix(rob.est1[,1], nrow = m)
            for(k in seq_len(m)){
                IDX <- INDEX[[k]]
                l.pm <- IDX[,1]
                if(ncol(IDX) == 2){
                    l.mm <- IDX[,2]
                }else{
                    l.mm <- integer()
                }
                pps.pm <- intensData[l.pm, , drop = FALSE]
                pps.mm <- intensData[l.mm, , drop = FALSE]
                pps.im <- pps.mm
                l <- t(t(pps.mm >= pps.pm) & (sb[k,] > contrast.tau))
                pps.im[l] <- t(t(pps.pm)/2^sb[k,])[l]
                l <- t(t(pps.mm >= pps.pm) & (sb[k,] <= contrast.tau))
                pps.im[l] <- t(t(pps.pm)/2^(contrast.tau/(1 + (contrast.tau - sb[k,])/scale.tau)))[l]
                pm.corrected <- matrix(pmax.int(pps.pm - pps.im, delta), ncol = ncol(pps.pm))
                colnames(pm.corrected) <- colnames(pps.pm)
                rownames(pm.corrected) <- rownames(pps.pm)
                res[[k]] <- pm.corrected
            }
            if(verbose) cat(" done.\n")
        }else{
            if(verbose) cat("Extract PM data ...")
            pm.only <- function(INDEX, x){
                l.pm <- INDEX[,1]
                x[l.pm, , drop = FALSE]
            }
            res <- lapply(INDEX, pm.only, x = intensData)
            if(verbose) cat(" done.\n")
        }
        if(verbose) cat("Computing expression values ...")
        for(k in nr){
            ind <- which(NROW == k)
            temp <- matrix(do.call(rbind, res[ind]), nrow = k)
            ind1 <-  as.vector(sapply(seq_len(n)-1, function(x, ind, m){ ind + x*m }, ind = ind, m = m))
            rob.est[ind1, 1:2] <- robloxbioc(log2(t(temp)), eps = eps, eps.lower = eps.lower, 
                                             eps.upper = eps.upper, steps = steps, fsCor = fsCor, 
                                             mad0 = mad0)
            rob.est[ind1, 2] <- rob.est[ind1, 2]/sqrt(k)
        }
        if(verbose) cat(" done.\n")
        exp.mat <- 2^matrix(rob.est[,1], nrow = m) + add.constant
        se.mat <- 2^matrix(rob.est[,2], nrow = m)

        dimnames(exp.mat) <- list(ids, sampleNames(x))
        dimnames(se.mat) <- list(ids, sampleNames(x))
        eset <- new("ExpressionSet", phenoData = phenoData(x), 
            experimentData = experimentData(x), exprs = exp.mat, 
            se.exprs = se.mat, annotation = annotation(x))
        if(normalize){
            if(verbose) cat("Scale normalization ...")
            for (i in 1:ncol(exprs(eset))) {
                slg <- exprs(eset)[, i]
                sf <- sc/mean(slg, trim = 0.02)
                reported.value <- sf * slg
                exprs(eset)[, i] <- reported.value
            }
            if(verbose) cat(" done.\n")
        }
        return(eset)
    })
