print.OBsProb <-
function (x, X = TRUE, resp = TRUE, factors = TRUE, models = TRUE, 
    nTop = 10, digits = 3, plt = FALSE, verbose = FALSE, Sh=TRUE, CV=TRUE, ...) 
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
                calc <- c(x$N, x$COLS, x$BLKS, x$MXFAC, x$MXINT, 
                x$mdcnt)
            names(calc) <- c("nRun", "nFac", "nBlk", "mFac", 
                "mInt", "totMod")
    out.list <- list(calc = calc)
    print(round(calc, digits = digits))
    if (plt) 
        plot.OBsProb(x, code = TRUE)
    if (factors) {
    cat("\n Factor probabilities:\n")
        prob <- data.frame(Factor = names(x$prob), 
            Prob = round(x$prob, digits), row.names = seq(length(x$prob)))
        print(prob, digits = digits)
        out.list[["probabilities"]] <- prob
    }
    if (models) {
        cat("\n Model probabilities:\n")
        ind <- seq(min(nTop, x$NTOP))
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
    if(Sh){
        cat("\n Shannon index:\n")
        Sh= -(sum(x$ptop[which(x$ptop>0)]*log(x$ptop[which(x$ptop>0)])))/log(length(x$ptop))
        print(Sh,digits=digits)
        out.list[["Sh"]] <- Sh  
    }
    if(CV){
        cat("\n CV:\n")
        CV= sd(x$prob[-1])*sqrt((x$MXFAC-1)/x$MXFAC)/mean(x$prob[-1])
        print(CV,digits=digits)
        out.list[["CV"]] <- CV  
    }
    invisible(out.list)
}