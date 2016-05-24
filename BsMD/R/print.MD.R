print.MD <-
function (x, X = FALSE, resp = FALSE, Xcand = TRUE, models = TRUE, 
    nMod = x$nMod, digits = 3, verbose = FALSE, ...) 
{
    if (verbose) {
        print(unclass(x))
        return(invisible(NULL))
    }
    nFac <- ncol(x$X) - x$blk
    if (X) {
        cat("\n Design Matrix:\n")
        print(x$X)
    }
    if (resp) {
        cat("\n Response vector:\n")
        cat(round(x$Y, digits = digits), fill = 80)
    }
    cat("\n Base:\n")
    calc <- c(x$N0, x$COLS, x$BL, x$CUT, x$GAMMA, x$GAM2, x$NM)
    names(calc) <- c("nRuns", "nFac", "nBlk", "maxInt", "gMain", 
        "gInter", "nMod")
    print(calc)
    cat("\n Follow up:\n")
    out <- c(x$N, x$NRUNS, x$ITMAX, x$NSTART)
    names(out) <- c("nCand", "nRuns", "maxIter", "nStart")
    print(out)
    calc <- c(calc, out)
    out.list <- list(calc = calc)
    if (models && x$NM > 0) {
        cat("\n Competing Models:\n")
        ind <- seq(x$NM)
        Prob <- round(x$P, digits)
        NumFac <- x$NF
        Sigma2 <- round(x$SIGMA2, digits)
        Factors <- apply(x$JFAC, 1, function(x) ifelse(all(x == 
            0), "none", paste(x[x != 0], collapse = ",")))
        dd <- data.frame(Prob, Sigma2, NumFac, Factors)
        print(dd, digits = digits, right = FALSE)
        out.list[["models"]] <- dd
    }
    if (Xcand) {
        cat("\n Candidate runs:\n")
        print(round(x$Xcand, digits))
    }
    if (any(x$D <= 0)) 
        ind <- min(which(x$D <= 0))
    else ind <- x$NTOP
    toprun <- data.frame(D = x$TOPD, x$TOPDES)
    ind <- min(nMod, ind)
    cat("\n   Top", ind, "runs:\n")
    print(dd <- round(toprun[seq(ind), ], digits))
    out.list[["follow.up"]] <- dd
    invisible(out.list)
}
