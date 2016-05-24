ordtest <- function (ord,var,dim=1:ncol(ord$points),index='euclidean',nitr=1000)
{
    if (inherits(ord, c("pco", "nmds", "metaMDS"))) {
        points <- ord$points
    } else if (inherits(ord, 'pca')) {
        points <- ord$scores
    } else {
        stop('ordtest is only defined for pca, pco, nmds, and metaMDS objects')
    }
    tdist <- 0
    observed <- 0
    reps <- rep(0,nitr-1)
    variable <- deparse(substitute(var))
    var <- factor(var)
    for (i in levels(var)) {
        mask <- var == i
        tdist <- tdist + sum(dist(points[mask,dim],index))
    }
    observed <- tdist
    for (i in 1:(nitr-1)) {
        tdist <- 0
        var <- sample(var,length(var),replace=FALSE)
        for (j in levels(var)) {
            mask <- var == j
            tdist <- tdist + sum(dist(points[mask,dim],index))
        }
        reps[i] <- tdist
    }
    pval <- (sum(reps<=observed)+1)/nitr
    print(paste(variable,'<',pval))
    out <- list(obs=observed,p=pval,reps=reps)
    invisible(out)
}
