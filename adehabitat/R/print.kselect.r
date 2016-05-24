"print.kselect" <- function (x, ...)
{
    cat("Duality diagramm\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$nf:", x$nf, "axis-components saved")
    cat("\n$rank: ")
    cat(x$rank)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n")
    else cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("vector", "length",
                                            "mode", "content")))
    sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[4, ] <- c("$initfac", length(x$initfac),
                    mode(x$initfac), "initial factor")
    sumry[5, ] <- c("$initwei", length(x$initwei),
                    mode(x$initwei), "row weights of inittab")
    class(sumry) <- "table"
    print(sumry, ...)
    cat("\n")
    sumry <- array("", c(10, 4), list(1:10, c("data.frame", "nrow",
                                            "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[6, ] <- c("$initab", nrow(x$initab),
                    ncol(x$initab), "initial table centered per animal")
    sumry[7, ] <- c("$as", nrow(x$as), ncol(x$as), "axis upon kselect axis")
    sumry[8, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "rows of initab upon kselect axis")
    sumry[9, ] <- c("$mus", nrow(x$mus), ncol(x$mav), "mean use on kselect axis")
    sumry[10, ] <- c("$mav", nrow(x$mus), ncol(x$mav), "mean available on kselect axis")
    class(sumry) <- "table"
    print(sumry)
}

