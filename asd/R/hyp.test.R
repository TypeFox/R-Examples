hyp.test <-
function (comb.test, level = level, full.hyp = FALSE) 
{
    crit <- qnorm(1 - level)
    treats <- length(comb.test$hyp.comb)
    zscores <- comb.test$zscores
    hyp.reject <- matrix(0, nrow = 1, ncol = treats)
    colnames(hyp.reject) <- colnames(comb.test$zscores[[1]])
    rownames(hyp.reject) <- 1
    hyp.reject <- matrix(0, nrow = 1, ncol = treats)
    colnames(hyp.reject) <- colnames(zscores[[1]])
    rownames(hyp.reject) <- 1
    vhyp.comb <- unlist(comb.test$hyp.comb)
    hyp.len <- rep(1:treats, times = choose(treats, 1:treats))
    vzscores <- rep(unlist(comb.test$zscores), hyp.len)
    matzscores <- matrix(NA, nrow = length(vzscores)/treats, 
        ncol = treats)
    matrejects <- matzscores
    for (i in 1:treats) {
        matzscores[, i] <- vzscores[vhyp.comb == i]
        if (full.hyp == TRUE) {
            matrejects[, i] <- as.numeric(matzscores[, i] >= 
                crit)
        }
        if (sum(as.numeric(matzscores[, i] >= crit)) == length(matzscores[, 
            i])) {
            hyp.reject[i] <- 1
        }
    }
    colnames(matzscores) <- colnames(comb.test$zscores[[1]])
    colnames(matrejects) <- colnames(comb.test$zscores[[1]])
    rownames(matzscores) <- 1:(length(vzscores)/treats)
    rownames(matrejects) <- 1:(length(vzscores)/treats)
    if (full.hyp == TRUE) {
        hyp.names <- NULL
        for (i in 1:treats) {
            hyp.names <- append(hyp.names, colnames(comb.test$hyp.comb[[i]]))
        }
        vhypcomb <- rep(hyp.names, hyp.len)
        mathypcomb <- matrix(NA, nrow = length(vzscores)/treats, 
            ncol = treats)
        for (i in 1:treats) {
            mathypcomb[, i] <- vhypcomb[vhyp.comb == i]
        }
        colnames(mathypcomb) <- colnames(comb.test$zscores[[1]])
        rownames(mathypcomb) <- 1:(length(vzscores)/treats)
    }
    if (full.hyp == TRUE) {
        list(reject = hyp.reject, all.rejects = matrejects, all.hyp = mathypcomb)
    }
    else {
        list(reject = hyp.reject)
    }
}
