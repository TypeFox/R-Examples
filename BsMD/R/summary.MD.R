summary.MD <-
function (object, digits = 3, verbose = FALSE, ...) 
{
    if (verbose) {
        print(unclass(object))
        return(invisible(NULL))
    }
    nFac <- ncol(object$X) - object$blk
    cat("\n Base:\n")
    calc <- c(object$N0, object$COLS, object$BL, object$CUT, 
        object$GAMMA, object$GAM2, object$NM)
    names(calc) <- c("nRuns", "nFac", "nBlk", "maxInt", "gMain", 
        "gInter", "nMod")
    print(calc)
    cat("\n Follow up:\n")
    out <- c(object$N, object$NRUNS, object$ITMAX, object$NSTART)
    names(out) <- c("nCand", "nRuns", "maxIter", "nStart")
    print(out)
    calc <- c(calc, out)
    out.list <- list(calc = calc)
    if (any(object$D <= 0)) 
        ind <- min(which(object$D <= 0))
    else ind <- object$NTOP
    toprun <- data.frame(D = object$TOPD, object$TOPDES)
    ind <- min(10, ind)
    cat("\n   Top", ind, "runs:\n")
    print(dd <- round(toprun[seq(ind), ], digits))
    out.list[["follow.up"]] <- dd
    invisible(out.list)
}
