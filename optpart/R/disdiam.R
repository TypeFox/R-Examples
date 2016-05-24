disdiam <- function(x,dist,digits)
{
    UseMethod("disdiam")
}

disdiam.default <- function (x, dist, digits = 3)
{
    clustering <- clustify(x)
    cname <- levels(clustering)
    csize <- length(levels(clustering))
    nclustering <- as.numeric(clustering)

    if (class(dist) != "dist")
        stop("The second argument must an object of class dist")
    dist <- as.matrix(dist)

    diam <- rep(0, csize)
    for (i in 1:csize) {
        pnt <- nclustering == i
        subdis <- dist[pnt, pnt]
        diam[i] <- max(subdis)
    }
    diam <- as.numeric(format(diam, digits = digits, nsmall = digits))
    res <- data.frame(cname, csize, diam)
    names(res) <- c("cluster", "N", "diameter")
    mean <- sum(res$N[res$N > 1] * res$diameter[res$N > 1])/sum(res$N[res$N > 1])
    out <- list(diameters = res, mean = mean)
    class(out) <- "disdiam"
    out
}

disdiam.stride <- function (x, dist, digits = 3)
{
    if (class(x) != 'stride') 
        stop('You must pass an object of class stride')
    res <- rep(NA, ncol(x$clustering))
    dist <- as.matrix(dist)
    for (i in 1:ncol(x$clustering)) {
        members <- table(x$clustering[,i])
        sum <- 0
        for (j in 1:x$seq[i]) {
            if (members[j] > 1) {
                pnt <- x$clustering[,i] == j
                subdis <- dist[pnt,pnt]
                diam <- max(subdis)
                sum <- sum + diam * members[j]
            }
        }
    res[i] <- round(sum/sum(members[members>1]),digits)
    }
    clusters <- x$seq
    diameters <- res
    out <- data.frame(clusters, diameters)
    out
}

print.disdiam <- function(x, ...)
{
    print(x$diameters)
    cat(paste('\nMean = ',format(x$mean,digits=4),"\n"))
}

