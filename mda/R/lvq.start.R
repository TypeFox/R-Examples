lvq.start <-
function (x, g, subclasses) 
{
    cnames <- levels(fg <- factor(g))
    J <- length(cnames)
    g <- as.numeric(g)
    weights <- as.list(cnames)
    names(weights) <- cnames
    subclasses <- rep(subclasses, length = length(cnames))
    size <- sum(subclasses)
    cb <- lvqinit(x, g, size = size)
    TT <- olvq1(x, g, codebk = cb)
    TT <- lvq3(x, g, codebk = TT)
    cl <- as.numeric(TT$cl)
    R <- length(cl)
    cx <- TT$x
    p <- ncol(cx)
    for (j in seq(J)) {
        which <- cl == j
        number <- sum(which)
        if (number == 0) {
            cx <- rbind(cx, apply(x[g == j, ], 2, mean))
            cl <- c(cl, j)
            wmj <- matrix(1, sum(g == j), 1)
            number <- 1
        }
        else if (number == 1) 
            wmj <- matrix(1, sum(g == j), 1)
        else {
            jcx <- cx[which, ]
            jcl <- seq(number)
            jcluster <- lvqtest(list(x = jcx, cl = jcl), x[g == 
                j, ])
            needed <- unique(jcluster)
            rcl <- rep(0, number)
            rcl[needed] <- j
            cl[which] <- rcl
            wmj <- diag(number)[jcluster, needed, drop = FALSE]
            number <- length(needed)
        }
        dimnames(wmj) <- list(NULL, paste("s", seq(number), sep = ""))
        weights[[j]] <- wmj
    }
    TT <- cl > 0
    list(x = cx[TT, , drop = FALSE], cl = factor(cl[TT], labels = cnames), 
        weights = weights)
}

