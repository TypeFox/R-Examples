RCconstrains <-
function (ncategories, homogeneous) 
{
    ncategories1 <- ncategories - 1
    nodf <- helpvec <- rep.int(0, ncategories1)
    for (i in 1:(ncategories1 - 1)) {
        helpmat <- t(combn(c(1:ncategories1), i))
        for (j in 1:nrow(helpmat)) {
            helpvec <- rep.int(0, ncategories1)
            helpvec[helpmat[j, ]] <- 1
            nodf <- rbind(nodf, helpvec)
        }
    }
    nodf <- unique(nodf)
    n1 <- nrow(nodf)
    parscores <- matrix(1:ncategories, nrow = n1, ncol = ncategories, 
        byrow = TRUE)
    for (j in 1:ncategories1) {
        parscores[nodf[, j] == 1, j + 1] <- parscores[nodf[, 
            j] == 1, j]
    }
    parscores <- t(apply(parscores, 1, function(x) as.numeric(factor(x))))
    ans <- list(parscores = parscores, nodf = as.numeric(rowSums(nodf)))
    if (!homogeneous) {
        Homogeneous <- ans
        n1 <- length(Homogeneous$nodf)
        parscores <- cbind(apply(Homogeneous$parscores, 2, function(x) rep.int(x, 
            n1)), apply(Homogeneous$parscores, 2, function(x) rep(x, 
            each = n1)))
        nodf <- rep(Homogeneous$nodf, each = n1) + rep.int(Homogeneous$nodf, 
            n1)
        orderedindices <- order(nodf)
        ans <- list(parscores = parscores[orderedindices, ], 
            nodf = nodf[orderedindices])
    }
    ans
}

