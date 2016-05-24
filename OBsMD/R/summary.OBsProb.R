summary.OBsProb <-
function (object, nTop = 10, digits = 3, ...) 
{
    nFac <- ncol(object$X) - object$blk
    cat("\n Calculations:\n")
            calc <- c(object$N, object$COLS, object$BLKS, object$MXFAC, 
                object$MXINT, object$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "totMod")
    out.list <- list(calc = calc)
    print(round(calc, digits = digits))
    prob <- data.frame(Factor = names(object$prob),  
        Prob = round(object$prob, digits), row.names = seq(length(object$prob)))
        cat("\n Factor probabilities:\n")
        print(prob, digits = digits)
        cat("\n Model probabilities:\n")
        ind <- seq(min(nTop, object$NTOP))
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
     invisible(out.list)
}
