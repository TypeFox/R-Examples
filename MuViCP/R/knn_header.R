#Functions of general utility for knn-type programs
eu.dist <- function(x,y)          sum((x-y)^2)

ab.dist <- function(x,y)          sum(abs(x-y))

mh.dist <- function(x, y, A)      t(x-y) %*% solve(A) %*% (x-y)

eu.dist.matY <- function(x, Y)    apply(Y, MARGIN = 1, FUN = function(y) eu.dist(x, y))

ab.dist.matY <- function(x, Y)    apply(Y, MARGIN = 1, FUN = function(y) ab.dist(x, y))

mh.dist.matY <- function(x, Y, A) apply(Y, MARGIN = 1, FUN = function(y) mh.dist(x, y, A))

least.k <- function(x, k)         x[order(x)[1:k]]

which.least.k <- function(x, k)     names(x)[(order(x))[1:k]]

least.p <- function(x, p = 0.05)       x[order(x)[1:floor(p * length(x))]]

which.least.p <- function(x, p = 0.05)   names(x)[order(x)[1:floor(p * length(x))]]

##P is projected data
get.NN <- function(P, k = 2, p = !k,
                   test, train,
                   dist.type = c('euclidean', 'absolute', 'mahal'),
                   nn.type   = c('which', 'dist', 'max'))
    {
        tmp <- switch(dist.type,
                      euclidean = apply(P[test,], MARGIN = 1, FUN = eu.dist.matY, Y = P[train,]),
                      absolute  = apply(P[test,], MARGIN = 1, FUN = ab.dist.matY, Y = P[train,]),
                      mahal     = apply(P[test,], MARGIN = 1, FUN = mh.dist.matY, Y = P[train,], A = cov(P)))
        row.names(tmp) <- train
        colnames(tmp) <- test
        if(p != FALSE) k <- floor(p * length(train))
        nears <- switch(nn.type,
                        which = matrix(as.numeric(apply(tmp, MARGIN = 2, which.least.k, k)), ncol = length(test)),
                        dist  = apply(tmp, MARGIN = 2, least.k, k),
                        max   = apply(tmp, MARGIN = 2, least.k, k)[k,])
        return(nears)
    }
