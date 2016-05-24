summary.BsProb <-
function (object, nMod = 10, digits = 3, ...) 
{
    nFac <- ncol(object$X) - object$blk
    cat("\n Calculations:\n")
    if (object$INDGAM == 0) {
        if (object$INDG2 == 0) {
            calc <- c(object$N, object$COLS, object$BLKS, object$MXFAC, 
                object$MXINT, object$P, object$GAMMA, object$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "p", "g", "totMod")
        }
        else {
            calc <- c(object$N, object$COLS, object$BLKS, object$MXFAC, 
                object$MXINT, object$P, object$GAMMA[1], object$GAMMA[2], 
                object$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "p", "g[main]", "g[int]", "totMod")
        }
    }
    else {
        calc <- c(object$N, object$COLS, object$BLKS, object$MXFAC, 
            object$MXINT, object$P, object$GAMMA[1], object$GAMMA[object$NGAM], 
            object$mdcnt)
        names(calc) <- c("nRun", "nFac", "nBlk", "mFac", "mInt", 
            "p", "g[1]", paste("g[", object$NGAM, "]", sep = ""), 
            "totMod")
    }
    out.list <- list(calc = calc)
    print(round(calc, digits = digits))
    prob <- data.frame(Factor = names(object$sprob), Code = rownames(object$prob), 
        Prob = round(object$sprob, digits), row.names = seq(length(object$sprob)))
    if (object$INDGAM == 0) {
        cat("\n Factor probabilities:\n")
        print(prob, digits = digits)
        cat("\n Model probabilities:\n")
        ind <- seq(min(nMod, object$NTOP))
        Prob <- round(object$ptop, digits)
        NumFac <- object$nftop
        Sigma2 <- round(object$sigtop, digits)
        Factors <- apply(object$jtop, 1, function(x) ifelse(all(x == 
            0), "none", paste(x[x != 0], collapse = ",")))
        dd <- data.frame(Prob, Sigma2, NumFac, Factors)[ind, 
            ]
        print(dd, digits = digits, right = FALSE)
        out.list[["probabilities"]] <- prob
        out.list[["models"]] <- dd
    }
    if (object$INDGAM == 1) {
        cat("\n Posterior probabilities for each gamma value:\n")
        print(dd <- round(rbind(gamma = object$GAMMA, object$prob), 
            digits = digits))
        out.list[["probabilities"]] <- dd
    }
    invisible(out.list)
}
