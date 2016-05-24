"print.enfa" <- function (x, ...)
{
    if (!inherits(x, "enfa"))
        stop("Object of class 'enfa' expected")
    cat("ENFA")
    cat("\n$call: ")
    print(x$call)
    cat("\nmarginality: ")
    cat(signif(x$m, 4))
    cat("\neigen values of specialization: ")
    l0 <- length(x$s)
    cat(signif(x$s, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...")
    cat("\n$nf:", x$nf, "axis of specialization saved")
    cat("\n")
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("vector", "length",
                                            "mode", "content")))
    sumry[1, ] <- c("$pr", length(x$pr), mode(x$pr), "vector of presence")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[4, ] <- c("$mar", length(x$mar),
                    mode(x$mar), "coordinates of the marginality vector")
    sumry[5, ] <- c("$s", length(x$s),
                    mode(x$s), "eigen values of specialization")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(3, 4), list(1:3, c("data.frame", "nrow",
                                            "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    class(sumry) <- "table"
    print(sumry)
    if (length(names(x)) > 11) {
        cat("\nother elements: ")
        cat(names(x)[12:(length(x))], "\n")
    }
}

