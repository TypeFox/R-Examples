print.BsProb <-
function (x, X = TRUE, resp = TRUE, factors = TRUE, models = TRUE, 
    nMod = 10, digits = 3, plt = FALSE, verbose = FALSE, ...) 
{
    if (verbose) {
        print(unclass(x))
        return(invisible(NULL))
    }
    nFac <- ncol(x$X) - x$blk
    if (X) {
        cat("\n Design Matrix:\n")
        print(round(x$X, digits))
    }
    if (resp) {
        cat("\n Response vector:\n")
        cat(round(x$Y, digits = digits), fill = 80)
    }
    cat("\n Calculations:\n")
    if (x$INDGAM == 0) {
        if (x$INDG2 == 0) {
            calc <- c(x$N, x$COLS, x$BLKS, x$MXFAC, x$MXINT, 
                x$P, x$GAMMA, x$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "p", "g", "totMod")
        }
        else {
            calc <- c(x$N, x$COLS, x$BLKS, x$MXFAC, x$MXINT, 
                x$P, x$GAMMA[1], x$GAMMA[2], x$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "p", "g[main]", "g[int]", "totMod")
        }
    }
    else {
        calc <- c(x$N, x$COLS, x$BLKS, x$MXFAC, x$MXINT, x$P, 
            x$GAMMA[1], x$GAMMA[x$NGAM], x$mdcnt)
        names(calc) <- c("nRun", "nFac", "nBlk", "mFac", "mInt", 
            "p", "g[1]", paste("g[", x$NGAM, "]", sep = ""), 
            "totMod")
    }
    out.list <- list(calc = calc)
    print(round(calc, digits = digits))
    if (plt) 
        plot.BsProb(x, code = TRUE)
    if (factors) {
        if (x$INDGAM == 1) 
            cat("\n Weighted factor probabilities:\n")
        else cat("\n Factor probabilities:\n")
        prob <- data.frame(Factor = names(x$sprob), Code = rownames(x$prob), 
            Prob = round(x$sprob, digits), row.names = seq(length(x$sprob)))
        print(prob, digits = digits)
        out.list[["probabilities"]] <- prob
    }
    if (x$INDGAM == 0 & models) {
        cat("\n Model probabilities:\n")
        ind <- seq(min(nMod, x$NTOP))
        Prob <- round(x$ptop, digits)
        NumFac <- x$nftop
        Sigma2 <- round(x$sigtop, digits)
        Factors <- apply(x$jtop, 1, function(x) ifelse(all(x == 
            0), "none", paste(x[x != 0], collapse = ",")))
        dd <- data.frame(Prob, Sigma2, NumFac, Factors)[ind, 
            ]
        print(dd, digits = digits, right = FALSE)
        out.list[["models"]] <- dd
    }
    if (x$INDGAM == 1) {
        cat("\n Values of posterior density of gamma:\n")
        dd <- data.frame(gamma = x$GAMMA, pgam = x$pgam)
        out.list[["gamma.density"]] <- dd
        print(dd, digits = digits)
        cat("\n Posterior probabilities for each gamma value:\n")
        print(dd <- round(rbind(gamma = x$GAMMA, x$prob), digits = digits))
        out.list[["probabilities"]] <- dd
    }
    invisible(out.list)
}
