canomi <- function(dudiX, Y, scannf = TRUE, nf = 2)
{
    if (!inherits(dudiX, "dudi"))
        stop("Object of class dudi expected")
    lig1 <- nrow(dudiX$tab)
    if (!is.data.frame(Y))
        stop("Y is not a data.frame")
    lig2 <- nrow(Y)
    if (lig1 != lig2)
        stop("Non equal row numbers")
    w1 <- apply(Y, 2, sum)
    if (any(w1 <= 0))
        stop(paste("Column sum <=0 in Y"))
    Z <- as.matrix(dudiX$tab)
    Y <- as.matrix(Y)
    lw <- dudiX$lw
    cw <- dudiX$cw
    lw <- lw/sum(lw)
    F <- apply(Y,2,function(x) x/sum(x))

    ## matrix of cross-product
    M <- crossprod(F,Z)

    ## Inverse of ZtRZ
    ZtRZ <- crossprod(apply(Z,2,function(x) x*sqrt(lw)))
    eS <- eigen(ZtRZ)
    ZtRZm12 <- eS$vectors[,eS$values>1e-7] %*%
        diag(1/sqrt(eS$values[eS$va>1e-7])) %*%
            t(eS$vectors[,eS$va>1e-7])
    Dt <- diag(apply(Y,2,sum)/sum(Y))
    C <- ZtRZm12%*%t(M)%*%Dt%*%M%*%ZtRZm12

    ## Calculation of the matrix to diagonalize
    res <- list()
    res$tab <- as.data.frame(M)
    eig1 <- eigen(C, symmetric = TRUE)
    eig <- eig1$values
    rank <- sum((eig/eig[1]) > 1e-7)
    if (scannf) {
        barplot(eig[1:rank])
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0)
        nf <- 2
    if (nf > rank)
        nf <- rank

    ## The eigenvalues
    res$eig <- eig[1:rank]

    ## The rank
    res$rank <- rank

    ## The number of axes kept
    res$nf <- nf

    dval <- sqrt(res$eig)[1:nf]
    auxi <- eig1$vectors[, 1:nf]
    auxi2 <- data.frame(M %*% ZtRZm12%*%auxi)

    ## The axes of the analyse
    auxi <- data.frame(auxi)
    names(auxi) <- paste("CS", (1:nf), sep = "")
    row.names(auxi) <- make.names(names(res$tab), unique = TRUE)
    res$c1 <- auxi

    ## The projection of the marginality vectors on the axes of the analysis
    names(auxi2) <- paste("Axis", (1:nf), sep = "")
    row.names(auxi2) <- row.names(res$tab)
    res$li <- auxi2

    ## The column metric
    res$cm <- ZtRZm12

    ## The projection of the RUs on the axes of the analysis
    res$ls <- Z%*%ZtRZm12%*%as.matrix(res$c1)

    ## The projection of the axes of the preliminary dudi.*
    ## on the axes of the analysis
    U <- as.matrix(res$c1)
    U <- data.frame(t(as.matrix(dudiX$c1)) %*% U)
    row.names(U) <- names(dudiX$li)
    names(U) <- names(res$li)
    res$as <- U

    ## The correlation between environmental variables and the
    ## scores of the RUs
    res$cor <- t(as.matrix(Z))%*%apply(as.matrix(res$ls), 2,function(x) x*lw)
    res$call <- match.call()
    class(res) <- "canomi"
    return(res)
}


print.canomi <- function (x, ...)
{
    if (!inherits(x, "canomi"))
        stop("to be used with 'canomi' object")
    cat("Canonical OMI analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$rank (rank)     :", x$rank)
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                          "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")

    sumry <- array("", c(7, 4), list(1:7, c("data.frame", "nrow",
                                            "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "crossed array (averaging animals/sites)")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "coordinates of the marginality vectors")
    sumry[3, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "variables normed scores")
    sumry[4, ] <- c("$cor", nrow(x$cor), ncol(x$cor), "variables correlations")
    sumry[5, ] <- c("$cm", nrow(x$cm), ncol(x$cm), "variables metric for the analysis")
    sumry[6, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "resource units coordinates")
    sumry[7, ] <- c("$as", nrow(x$as), ncol(x$as), "axis upon niche axis")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}


plot.canomi <- function (x, xax = 1, yax = 2, ...)
{
    if (!inherits(x, "canomi"))
        stop("Use only with 'canomi' objects")
    if (x$nf == 1) {
        warnings("One axis only : not yet implemented")
        return(invisible())
    }
    if (xax > x$nf)
        stop("Non convenient xax")
    if (yax > x$nf)
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), respect = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.corcircle(x$as, xax, yax, sub = "Axis", csub = 2, clabel = 1.25)
    s.arrow(x$c1, xax, yax, sub = "Variable Scores", csub = 2, clabel = 1.25)
    s.arrow(x$cor, xax, yax, sub = "Correlations", csub = 2, clabel = 1.25)
    s.label(x$ls, xax, yax, clabel = 0, cpoint = 2, sub = "RUs and animals",
        csub = 2)
    s.distri(x$ls, eval(as.list(x$call)[[3]], sys.frame(0)), add.plot=TRUE)
    s.label(x$li, xax, yax, clabel = 1.5, add.plot = TRUE)
    s.arrow(x$li, sub = "Marginality vectors",
            csub = 2)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))

}
