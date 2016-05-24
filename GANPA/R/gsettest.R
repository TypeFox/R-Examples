GSE.Test.Main <-
function (gExprs.obj, gsets, gNET, 
    check.exprs = TRUE, msp.groups, size.min = 15, 
    size.max = 500, permN = 1000, randN = 30, permFDR.cutoff = 0.5, 
    output.label = "dataset", msp.correction = TRUE) 
{
    gsets.w.orig <- weight.gsets.test(gsets = gsets, isets = gNET)
    if (msp.correction) {
        gNET.multi <- create.sets.multi(sets = gNET, msp.groups = msp.groups)
        gsets.w.multi <- weight.gsets.with.msprot(gsets = gsets, 
            isets.multi = gNET.multi, msp.groups = msp.groups)
    }
    T.obj <- obtainT(gExprs = gExprs.obj, reverse = FALSE, check.exprs = check.exprs, 
        permN = permN)
    gsets.fil <- preproc.gsets(gsets = gsets, gbg = rownames(gExprs.obj$gExprs), 
        size.max = size.max, size.min = size.min)
    gsets.w.orig.fil <- gsets.w.orig[names(gsets.fil)]
    GSET.MeanAbs.NoW <- GSE.Test(T.obj = T.obj, gsets.weight = gsets.w.orig.fil, 
        method = "absmean", randN = randN)
    GSET.MeanAbs.OrigW <- GSE.Test(T.obj = T.obj, gsets.weight = gsets.w.orig.fil, 
        method = "abswmean", randN = randN)
    if (msp.correction) {
        gsets.w.multi.fil <- gsets.w.multi[names(gsets.fil)]
        GSET.MeanAbs.MultiW <- GSE.Test(T.obj = T.obj, gsets.weight = gsets.w.multi.fil, 
            method = "abswmean", randN = randN)
    }
    write.csv(summary.GSE.Test(GSET.MeanAbs.NoW, fdr.perm.cut = permFDR.cutoff), 
        file = paste(output.label, "MeanAbs.NoW.csv", sep = "."))
    write.csv(summary.GSE.Test(GSET.MeanAbs.OrigW, fdr.perm.cut = permFDR.cutoff), 
        file = paste(output.label, "MeanAbs.OrigW.csv", sep = "."))
    if (msp.correction) 
        write.csv(summary.GSE.Test(GSET.MeanAbs.MultiW, fdr.perm.cut = permFDR.cutoff), 
            file = paste(output.label, "MeanAbs.MultiW.csv", 
                sep = "."))
}
GSE.Test <-
function (T.obj, gsets.weight, method = "abswmean", randN = 1000, 
    restand = "perm") 
{
    weighted <- FALSE
    if (is.null(gs.names <- names(gsets.weight))) 
        names(gsets.weight) <- gs.names <- paste("gs", 1:length(gsets.weight), 
            sep = "_")
    if (!is.numeric(gsets.weight[[1]])) 
        stop("'gsets.weight' should be named numeric vectors.")
    if (method == "wmean" | method == "abswmean") 
        weighted <- TRUE
    if (method == "absmean" | method == "abswmean") {
        alternative <- "greater"
    }
    else if (method == "mean" | method == "wmean") 
        alternative <- "two.sided"
    tG.gsets <- unique(unlist(lapply(gsets.weight, names)))
    tG.array <- names(T.obj$T.observed)
    if (length(tG.com <- intersect(tG.gsets, tG.array)) < length(tG.gsets) * 
        0.1) 
        stop("<10% of the gene set genes can be found in the microarray.")
    weight.gsets.index <- lapply(gsets.weight, function(w) {
        g.ind.array <- match(names(w), tG.array)
        ind <- which(!is.na(g.ind.array))
        names(w) <- g.ind.array
        w[ind]
    })
    gs.scores.obs <- score.func(weight.gsets.index = weight.gsets.index, 
        gT = T.obj$T.observed, method = method)
    print("randomization...")
    gs.scores.rand <- NULL
    Null <- lapply(1:randN, function(i) {
        gT.rand <- sample(T.obj$T.observed)
        scores.rand <- score.func(weight.gsets.index = weight.gsets.index, 
            gT = gT.rand, method = method)
        gs.scores.rand <<- cbind(gs.scores.rand, scores.rand)
        return(NULL)
    })
    p.rand <- apply(gs.scores.rand - gs.scores.obs, 1, function(x) {
        length(which(x >= 0))/randN
    })
    p.rand <- switch(alternative, greater = p.rand, two.sided = min(p.rand, 
        1 - p.rand))
    print("permutation...")
    gs.scores.perm <- apply(T.obj$T.permuted, 2, function(gT.perm) {
        score.func(weight.gsets.index = weight.gsets.index, gT = gT.perm, 
            method = method)
    })
    p.perm <- apply(gs.scores.perm - gs.scores.obs, 1, function(x) {
        length(which(x >= 0))/ncol(gs.scores.perm)
    })
    p.perm <- switch(alternative, greater = p.perm, two.sided = pmin(p.perm, 
        1 - p.perm))
    if (restand == "rand") {
        print("restandardization based on row randomization...")
        ref.median <- apply(gs.scores.rand, 1, median)
        ref.sd <- apply(gs.scores.rand, 1, sd)
    }
    else if (restand == "perm") {
        print("restandardization based on sample permutation...")
        ref.median <- apply(gs.scores.perm, 1, median)
        ref.sd <- apply(gs.scores.perm, 1, sd)
    }
    gs.scores.obs.norm <- (gs.scores.obs - ref.median)/ref.sd
    print("Done...")
    stat0 <- cbind(gs.size = unlist(lapply(weight.gsets.index, 
        length)), gs.scores = gs.scores.obs, gs.scores.norm = gs.scores.obs.norm, 
        reference.median = ref.median, reference.sd = ref.sd, 
        p.rand = p.rand, fdr.rand = fdr.est(p.rand), p.perm = p.perm, 
        fdr.perm = fdr.est(p.perm))
    list(stat0 = stat0, method = method, alternative = alternative, 
        gsets.weight = gsets.weight, restand = restand)
}
obtainT <-
function (gExprs, permN = 1000, log2ratio = FALSE, reverse = FALSE, 
    check.exprs = FALSE) 
{
    if (check.exprs) 
        gExprs <- checkExprs(gExprs)
    if (is.factor(gExprs$sampleinfo[, "status"])) 
        gExprs$sampleinfo[, "status"] <- as.character(gExprs$sampleinfo[, 
            "status"])
    status <- gExprs$sampleinfo[, "status"]
    if (is.null(names(status))) 
        names(status) <- rownames(gExprs$sampleinfo)
    status <- status[which(!is.na(status) & status != "")]
    sampleNames <- colnames(gExprs$gExprs)
    if (length(status.uni <- unique(status)) != 2) 
        stop("obtainT(): status number is not 2!")
    status.uni <- sort(status.uni)
    if (reverse == TRUE) 
        status.uni <- rev(status.uni)
    status1 <- status.uni[1]
    status2 <- status.uni[2]
    sample1 <- intersect(names(status[which(status == status1)]), 
        sampleNames)
    sample2 <- intersect(names(status[which(status == status2)]), 
        sampleNames)
    if (length(sample1) < 3 | length(sample2) < 3) 
        stop("obtainT(): either status group 1 or 2 has < 3 samples!")
    if ((length(sample1) + length(sample2)) < length(sampleNames) * 
        0.5) 
        warning("obtainT(): poor overlap between colnames of expression matrix and status names, please have a check.")
    Mat <- gExprs$gExprs[, c(sample1, sample2)]
    status <- c(rep(1, length(sample1)), rep(2, length(sample2)))
    calT <- function(Mat, status, log2ratio = FALSE) {
        m <- length(ind1 <- which(status == 1))
        n <- length(ind2 <- which(status == 2))
        mat1 <- Mat[, ind1]
        mat2 <- Mat[, ind2]
        rmean1 <- rowMeans(mat1)
        rmean2 <- rowMeans(mat2)
        if (log2ratio) 
            return(rmean2 - rmean1)
        ss1 <- rowSums((mat1 - rmean1)^2)
        ss2 <- rowSums((mat2 - rmean2)^2)
        tt <- (m + n - 2)^0.5 * (rmean2 - rmean1)/((1/m + 1/n) * 
            (ss1 + ss2))^0.5
        return(list(T = tt, df = m + n - 2))
    }
    if (log2ratio) {
        ratio <- calT(Mat = Mat, status = status, log2ratio = TRUE)
        return(list(T.observed = ratio, df = NULL, T.permuted = NULL))
    }
    if (is.null(permN)) {
        permN.potent <- factorial(length(sample1) + length(sample2))/(factorial(length(sample1)) * 
            factorial(length(sample2)))
        print(paste("Number of possible different permnutations: ", 
            permN.potent, sep = ""))
        if (permN.potent < 3000) 
            permN <- 1000
        else permN <- 3000
        print(paste("Running automatically chosen permutation number: ", 
            permN, sep = ""))
    }
    else print(paste("Running the required permutation number: ", 
        permN, sep = ""))
    status.names <- names(status)
    status.permuted <- lapply(1:permN, function(x) {
        status.perm <- sample(status)
        names(status.perm) <- status.names
        return(status.perm)
    })
    T.stat <- calT(Mat = Mat, status = status)
    names(T.stat) <- c("T.observed", "df")
    T.permuted <- lapply(1:length(status.permuted), function(i) {
        if (i%%100 == 0) 
            cat(paste("\tperms ", i, "\n", sep = ""))
        status <- status.permuted[[i]]
        T.obj <- calT(Mat = Mat, status = status)
        T.obj$T
    })
    rowNames <- names(T.permuted[[1]])
    T.permuted <- matrix(unlist(T.permuted), nrow = length(rowNames))
    rownames(T.permuted) <- rowNames
    group <- c(controlgroup = status1, testgroup = status2)
    return(c(T.stat, group = list(group), permN = permN, list(T.permuted = T.permuted)))
}
score.func <-
function (weight.gsets.index, gT, method = "abswmean") 
{
    meth <- c("mean", "absmean", "wmean", "abswmean")
    if (length(intersect(method, meth)) != 1) 
        stop(paste("'method' should be only one of the following: ", 
            paste(meth, collapse = ", "), sep = ""))
    if (method == "mean" | method == "absmean") 
        weight.gsets.index <- lapply(weight.gsets.index, function(w) {
            w0 <- rep(1, length(w))
            names(w0) <- names(w)
            return(w0)
        })
    if (method == "absmean" | method == "abswmean") 
        gT <- abs(gT)
    out <- unlist(lapply(weight.gsets.index, function(weight) {
        gset.index <- as.integer(names(weight))
        s <- sum(weight * gT[gset.index])/sum(weight)
    }))
    return(out)
}
summary.GSE.Test <-
function (gse.obj, fdr.perm.cut = 0.1, fdr.rand.cut = NULL, topN = NULL, 
    rank = "prs") 
{
    if (length(grep("^restand$", names(gse.obj))) == 1) {
        stat0 <- gse.obj$stat0[order(gse.obj$stat0[, "p.perm"], 
            -abs(gse.obj$stat0[, "gs.scores.norm"])), ]
    }
    else {
        if (rank == "prs") 
            stat0 <- gse.obj$stat0[order(gse.obj$stat0[, "p.perm"], 
                gse.obj$stat0[, "p.rand"], 1/abs(gse.obj$stat0[, 
                  "gs.scores"])), ]
        else if (rank == "ps") 
            stat0 <- gse.obj$stat0[order(gse.obj$stat0[, "p.perm"], 
                -abs(gse.obj$stat0[, "gs.scores"])), ]
        else stop("unknown 'rank' option!")
    }
    stat <- stat0
    if (is.null(topN)) {
        if (!is.null(fdr.perm.cut)) 
            stat <- stat[which(stat[, "fdr.perm"] <= fdr.perm.cut), 
                ]
        if (!is.null(fdr.rand.cut)) 
            stat <- stat[which(stat[, "fdr.rand"] <= fdr.rand.cut), 
                ]
    }
    else stat <- stat[1:topN, ]
    if (ncol(stat) == 6) 
        colnames(stat) <- c("Size", "Score", "randP", "randFDR", 
            "permP", "permFDR")
    else colnames(stat) <- c("Size", "S", "NS", "refMedian", 
        "refSd", "randP", "randFDR", "permP", "permFDR")
    round(stat, 3)
}
write.gct <-
function (gExprs, file.gct, file.cls) 
{
    status <- gExprs$sampleinfo[, "status"]
    if (is.null(names(status))) {
        names(status) <- sampleid <- rownames(gExprs$sampleinfo)
    }
    else sampleid <- names(status)
    ind.fil <- which(!is.na(hit <- match(colnames(gExprs$gExprs), 
        sampleid)))
    mat <- gExprs$gExprs[, ind.fil]
    colnames(mat) <- sampleid[hit[ind.fil]]
    status <- status[hit[ind.fil]]
    mat <- collapseProbes(Mat = mat)
    cat("#1.2\n", file = file.gct)
    cat(paste(dim(mat), collapse = "\t"), file = file.gct, append = TRUE)
    cat("\n", file = file.gct, append = TRUE)
    mat <- cbind(NAME = row.names(mat), DESCRIPTION = rep("na", 
        nrow(mat)), mat)
    write.table(mat, file = file.gct, sep = "\t", col.names = TRUE, 
        row.names = FALSE, quote = FALSE, append = TRUE)
    status.uni <- sort(unique(status))
    status.cls <- rep(NA, length(status))
    for (i in 1:length(status.uni)) status.cls[status == status.uni[i]] <- i
    cat(paste(length(status), length(status.uni), 1, sep = " "), 
        file = file.cls)
    cat("\n", file = file.cls, append = TRUE)
    cat(paste("#", paste(status.uni, collapse = " "), "\n", sep = ""), 
        file = file.cls, append = TRUE)
    cat(paste(status.cls, collapse = " "), file = file.cls, append = TRUE)
    cat("\n", file = file.cls, append = TRUE)
}
fdr.est <-
function (p) 
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}
