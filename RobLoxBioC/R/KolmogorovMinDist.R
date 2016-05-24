###############################################################################
## Methods for KolmogorovMinDist: computation of minimum Kolmogorov distance
###############################################################################
setMethod("KolmogorovMinDist", signature(x = "matrix",
                                         D = "Norm"),
    function(x, D, mad0 = 1e-4){
        stopifnot(is.numeric(x))
        ksdist <- function(pars, x){
            if(any(ina <- is.na(x))) x <- x[!ina]
            y <- .Internal(sort(x, FALSE))
            p.y <- pnorm(y, mean = pars[1], sd = pars[2])
            n <- length(y)
            max(p.y - (seq_len(n)-1)/n, seq_len(n)/n - p.y)
        }
        Med <- rowMedians(x, na.rm = TRUE)
        Mad <- pmax(rowMedians(abs(x-Med), na.rm = TRUE)/qnorm(0.75), mad0)
        startPars <- cbind(Med, Mad)
        M <- nrow(x)
        n <- ncol(x)
        res <- matrix(NA, nrow = M, ncol = 2)
        for(i in seq_len(M)){
            temp <- try(optim(par = startPars[i,], fn = ksdist, x = x[i,], 
                              method = "L-BFGS-B", lower = c(-Inf, 1e-15), 
                              upper = c(Inf, Inf))$value, silent = TRUE)
            if(!inherits(temp, "try-error")){
                res[i,1] <- temp
                res[i,2] <- n
            }
        }
        list("dist" = res[,1], "n" = res[,2])
    })

setMethod("KolmogorovMinDist", signature(x = "AffyBatch",
                                         D = "AbscontDistribution"),
    function(x, D, bg.correct = TRUE, pmcorrect = TRUE, verbose = TRUE){
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
        kmd <- numeric(m*n)
        ns <- numeric(m*n)
        if(pmcorrect){
            if(verbose) cat("Calculating PM/MM ...")
            div <- function(INDEX, x){
                l.pm <- INDEX[,1]
                if(ncol(INDEX) == 2)
                    l.mm <- INDEX[,2]
                else
                    l.mm <- integer()
                x[l.pm, , drop = FALSE]/x[l.mm, , drop = FALSE]
            }
            res <- lapply(INDEX, div, x = intensData)
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
        if(verbose) cat("Computing minimum Kolmogorov distance ...")
        for(k in nr){
            ind <- which(NROW == k)
            temp <- matrix(do.call(rbind, res[ind]), nrow = k)
            ind1 <-  as.vector(sapply(seq_len(n)-1, function(x, ind, m){ ind + x*m }, ind = ind, m = m))
            temp.res <- KolmogorovMinDist(log2(t(temp)), D = D)
            kmd[ind1] <- temp.res[["dist"]]
            ns[ind1] <- temp.res[["n"]]
        }
        if(verbose) cat(" done.\n")
        kmd.mat <- matrix(kmd, nrow = m)
        ns.mat <- matrix(ns, nrow = m)
        dimnames(kmd.mat) <- list(ids, sampleNames(x))
        dimnames(ns.mat) <- list(ids, sampleNames(x))
        list("dist" = kmd.mat, "n" = ns.mat)
    })

setMethod("KolmogorovMinDist", signature(x = "BeadLevelList",
                                         D = "AbscontDistribution"),
    function(x, D, log = FALSE, imagesPerArray = 1, what = "G", probes = NULL, arrays = NULL){
        BLData <- x
        arraynms <- arrayNames(BLData)
        if(!is.null(arrays) && !is.character(arrays)) arraynms <- arraynms[arrays]
        if(is.character(arrays)) arraynms <- which(arraynms %in% arrays)
        len <- length(arraynms)
        what <- match.arg(what, c("G", "R", "RG", "M", "A", "beta"))
        whatelse <- ""
        if(what == "RG"){
            if(BLData@arrayInfo$channels == "two"){
                what <- "G"
                whatelse <- "R"
            }else{
                stop("Need two-channel data to calculate summary R and G values")
            }
        }
        if(imagesPerArray == 1){
            sel <- getArrayData(BLData, what = "ProbeID", array = arraynms[1]) != 0
            pr <- getArrayData(BLData, what = "ProbeID", array = arraynms[1])[sel]
            finten <- getArrayData(BLData, what = what, log = log, array = arraynms[1])[sel]
            nasinf <- !is.finite(finten) | is.na(finten)
            finten <- finten[!nasinf]
        }
        else if(imagesPerArray == 2){
            if(length(arraynms)%%2 != 0) 
                stop("Need an even number of arrays when 'imagesPerArray=2'")
            arrayord <- order(arraynms)
            arraynms <- arraynms[arrayord]
            tmp <- unlist(strsplit(arraynms, "_"))
            chipnums <- tmp[seq(1, length(tmp), by = 3)]
            pos <- tmp[seq(2, length(tmp), by = 3)]
            stripnum <- as.numeric(tmp[seq(3, length(tmp), by = 3)])
            check <- ((chipnums[seq(1, len, by = 2)] == chipnums[seq(2, len, by = 2)]) 
                      & (pos[seq(1, len, by = 2)] == pos[seq(2, len, by = 2)]) 
                      & (stripnum[seq(1, len, by = 2)] == stripnum[seq(2, len, by = 2)] - 1))
            if(sum(check) != length(check)) stop("Missing arrays")
            sel1 <- getArrayData(BLData, what = "ProbeID", array = arraynms[1]) != 0
            sel2 <- getArrayData(BLData, what = "ProbeID", array = arraynms[2]) != 0
            pr <- append(getArrayData(BLData, what = "ProbeID", array = arraynms[1])[sel1], 
                        getArrayData(BLData, what = "ProbeID", array = arraynms[2])[sel2])
            finten <- append(getArrayData(BLData, what = what, log = log, array = arraynms[1])[sel1], 
                            getArrayData(BLData, what = what, log = log, array = arraynms[2])[sel2])
            nasinf <- !is.finite(finten) | is.na(finten)
            finten <- finten[!nasinf]
        }else{
            stop("You can only specify 1 or 2 images per array")
        }
        if(is.null(probes)) probes <- sort(unique(pr))
        probes <- probes[probes > 0 & !is.na(probes)]
        noprobes <- length(probes)
        pr <- pr[!nasinf]
        if (imagesPerArray == 1) {
            G.kmd <- matrix(0, nrow = noprobes, ncol = len)
            G.ns <- matrix(0, nrow = noprobes, ncol = len)
            colnames(G.kmd) <- colnames(G.ns) <- arraynms
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse == "R") 
                R.kmd <- R.ns <- G.kmd
            else R.kmd <- R.ns <- NULL
        }
        else if (imagesPerArray == 2) {
            G.kmd <- matrix(0, nrow = noprobes, ncol = (len/2))
            G.ns <- matrix(0, nrow = noprobes, ncol = (len/2))
            colnames(G.kmd) <- colnames(G.ns) <- arraynms[seq(1, len, by = 2)]
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse == "R") 
                R.kmd <- R.ns <- G.kmd
            else R.kmd <- R.ns <- NULL
        }
        i <- j <- 1
        while (j <= len) {
            probeIDs <- as.integer(pr)
            start <- 0
            G.res <- kmdBeadLevel(x = finten, D = D, probeIDs = probeIDs, probes = probes)
            G.kmd[,i] <- G.res[["dist"]]
            G.ns[,i] <- G.res[["n"]]
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[i]]]$R) && whatelse == "R") {
                if (imagesPerArray == 1) {
                    finten <- getArrayData(BLData, what = whatelse, log = log, array = arraynms[i])[sel]
                    nasinf <- !is.finite(finten) | is.na(finten)
                    finten <- finten[!nasinf]
                    binten <- rep(0, length(finten))
                }
                else if (imagesPerArray == 2) {
                    finten <- append(getArrayData(BLData, what = whatelse, log = log, array = arraynms[j])[sel1], 
                                    getArrayData(BLData, what = whatelse, log = log, array = arraynms[j + 1])[sel2])
                    nasinf <- !is.finite(finten) | is.na(finten)
                    finten <- finten[!nasinf]
                    binten <- rep(0, length(finten))
                }
                start <- 0
                R.res[,i] <- kmdBeadLevel(x = finten, D = D, probeIDs = probeIDs, probes = probes)
                R.kmd[,i] <- R.res[["dist"]]
                R.ns[,i] <- R.res[["n"]]
            }
            j <- j + imagesPerArray
            i <- i + 1
            rm(probeIDs)
            gc()
            if ((imagesPerArray == 1) && (i <= len)) {
                sel <- getArrayData(BLData, what = "ProbeID", array = arraynms[i]) != 0
                pr <- getArrayData(BLData, what = "ProbeID", array = arraynms[i])[sel]
                finten <- getArrayData(BLData, what = what, log = log, array = arraynms[i])[sel]
                nasinf <- !is.finite(finten) | is.na(finten)
                pr <- pr[!nasinf]
                finten <- finten[!nasinf]
            }
            else if ((imagesPerArray == 2) && (j < len)) {
                sel1 <- getArrayData(BLData, what = "ProbeID", array = arraynms[j]) != 0
                sel2 <- getArrayData(BLData, what = "ProbeID", array = arraynms[j + 1]) != 0
                pr <- append(getArrayData(BLData, what = "ProbeID", array = arraynms[j])[sel1], 
                             getArrayData(BLData, what = "ProbeID", array = arraynms[j + 1])[sel2])
                finten <- append(getArrayData(BLData, what = what, log = log, array = arraynms[j])[sel1], 
                                 getArrayData(BLData, what = what, log = log, array = arraynms[j + 1])[sel2])
                nasinf <- !is.finite(finten) | is.na(finten)
                pr <- pr[!nasinf]
                finten <- finten[!nasinf]
            }
        }
        if (whatelse == "R") {
            rownames(G.kmd) <- rownames(G.ns) <- rownames(R.kmd) <- rownames(R.ns) <- probes
            res <- list(G = list("dist" = G.kmd, "n" = G.ns), 
                        R = list("dist" = R.kmd, "n" = R.ns))
        }
        else {
            rownames(G.kmd) <- rownames(G.ns) <- probes
            res <- list("dist" = G.kmd, "n" = G.ns)
        }
        return(res)
    })
kmdBeadLevel <- function(x, D, probeIDs, probes){
    comIDs <- intersect(probeIDs, probes)
    x <- x[probeIDs %in% comIDs]
    probeIDs <- probeIDs[probeIDs %in% comIDs]
    noBeads <- as.vector(table(probeIDs))
    noBeads.uni <- as.integer(names(table(noBeads)))
    probes1 <- comIDs
    len1 <- length(probes1)
    kmd1 <- numeric(len1)
    ns1 <- numeric(len1)
    for(i in seq(along = noBeads.uni)){
        index <- noBeads == noBeads.uni[i]
        IDs <- probes1[index]
        if(noBeads.uni[i] == 1){
            kmd1[index] <- 0.5
        }else{
            temp <- matrix(x[probeIDs %in% IDs], ncol = noBeads.uni[i], byrow = TRUE)
            res <- KolmogorovMinDist(temp, D = D)
            kmd1[index] <- res[["dist"]]
            ns1[index] <- res[["n"]]
        }
    }    
    len <- length(probes)
    kmd <- numeric(len)
    ns <- numeric(len)
    nas <- !(probes %in% comIDs)
    kmd[nas] <- NA
    kmd[!nas] <- kmd1
    ns[nas] <- NA
    ns[!nas] <- ns1

    list("dist" = kmd, "n" = ns)
}
