print.rocJM <-
function (x, ...) {
    cat("\nAreas under the time-dependent ROC curves\n\n")
    cat("Estimation: Monte Carlo (", x$M, " samples)\n", sep = "")
    if (x$diffType == "absolute") {
        lx <- length(x$abs.diff)
        ld <- paste(round(x$abs.diff, 2), collapse = ", ")
    } else {
        lx <- length(x$rel.diff)
        ld <- paste(round(x$rel.diff, 2), collapse = ", ")
    }
    cat("Difference: ", x$diffType, ", lag = ", lx, 
        " (", ld,")\n", sep = "")
    cat("Thresholds range: (", round(x$min.cc, 2), ", ", 
        round(x$max.cc, 2), ")\n\n", sep = "")
    times <- x$times
    aucs <- x$AUCs
    for (i in seq_along(times)) {
        cat("Case:", names(times)[i], "\n")
        cat("Recorded time(s):", paste(round(x$times[[i]], 2), 
            collapse = ", "), "\n")
        ac <- if (is.matrix(aucs)) round(aucs[, i], 4) else 
            round(aucs[[i]], 4)
        thr <- round(x$optThr[[i]], 4)
        m <- cbind(x$dt, round(x$dt + 
            tail(x$times[[i]], 1), 2), ac, thr)
        colnames(m) <- if ((nc <- ncol(thr)) == 1) {
            c("dt", "t + dt", "AUC", "Cut")
        } else { 
            c("dt", "t + dt", "AUC", 
                paste("Cut.", 1:nc, sep = ""))
        }
        rownames(m) <- rep("", nrow(m))
        print(m)
        cat("\n")
    }
    invisible(x)
}
