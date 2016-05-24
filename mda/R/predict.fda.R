predict.fda <-
function (object, newdata, type = c("class", "variates", "posterior", 
    "hierarchical", "distances"), prior, dimension = J - 1, ...)
{
    dist <- function(x, mean, m = ncol(mean)) (scale(x, mean, 
        FALSE)^2) %*% rep(1, m)
    type <- match.arg(type)
    means <- object$means
    Jk <- dim(means)
    J <- Jk[1]
    k <- Jk[2]
    if (type == "hierarchical") {
        if (missing(dimension)) 
            dimension.set <- seq(k)
        else {
            dimension.set <- dimension[dimension <= k]
            if (!length(dimension.set)) 
                dimension.set <- k
            dimension <- max(dimension.set)
        }
    }
    else dimension <- min(max(dimension), k)
    if (missing(newdata)) 
        y <- predict(object$fit)
    else {
        if (inherits(newdata, "data.frame") || is.list(newdata)) {
            Terms <- delete.response(terms(object))
#            attr(Terms, "intercept") <- 0
            newdata <- model.matrix(Terms, newdata)
            if(attr(Terms,"intercept"))newdata=newdata[,-1,drop=FALSE]
          }
        y <- predict(object$fit, newdata)
    }
    y <- y %*% object$theta[, seq(dimension), drop = FALSE]
    lambda <- object$values
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    newdata <- scale(y, FALSE, sqima * alpha)
    if (missing(prior)) 
        prior <- object$prior
    else {
        if (any(prior < 0) | round(sum(prior), 5) != 1) 
            stop("innappropriate prior")
    }
    means <- means[, seq(dimension), drop = FALSE]
    switch(type, variates = return(newdata), class = {
        n <- nrow(newdata)
        prior <- 2 * log(prior)
        mindist <- dist(newdata, means[1, ], dimension) - prior[1]
        pclass <- rep(1, n)
        for (i in seq(2, J)) {
            ndist <- dist(newdata, means[i, ], dimension) - prior[i]
            l <- ndist < mindist
            pclass[l] <- i
            mindist[l] <- ndist[l]
        }
        ## 2001-10-27: Need to provide levels or else if we get an error
        ## if the predicted classes do no contain all possible classes.
        ## Reported by Greg Jefferis <jefferis@stanford.edu>, fix by
        ## Bj/orn-Helge Mevik <bjorn-helge.mevik@matforsk.no>.
        return(factor(pclass, levels = seq(J),
                      labels = dimnames(means)[[1]]))
    }, posterior = {
        pclass <- matrix(0, nrow(newdata), J)
        for (i in seq(J)) pclass[, i] <- exp(-0.5 * dist(newdata, means[i, 
            ], dimension)) * prior[i]
        dimnames(pclass) <- list(dimnames(newdata)[[1]], dimnames(means)[[1]])
        return(pclass/drop(pclass %*% rep(1, J)))
    }, hierarchical = {
        prior <- 2 * log(prior)
        Pclass <- vector("list", length(dimension.set))
        names(Pclass) <- paste("D", dimension.set, sep = "")
        for (ad in seq(along = dimension.set)) {
            d <- dimension.set[ad]
            dd <- seq(d)
            mindist <- dist(newdata[, dd, drop = FALSE], means[1, dd, drop = FALSE], 
                d) - prior[1]
            pclass <- rep(1, nrow(newdata))
            for (i in seq(2, J)) {
                ndist <- dist(newdata[, dd, drop = FALSE], means[i, dd, 
                  drop = FALSE], d) - prior[i]
                l <- ndist < mindist
                pclass[l] <- i
                mindist[l] <- ndist[l]
            }
            levels(pclass) <- dimnames(means)[[1]]
            Pclass[[ad]] <- pclass
        }
        rownames <- dimnames(newdata)[[1]]
        if (is.null(rownames)) 
            rownames <- paste(seq(nrow(newdata)))
        return(structure(Pclass, class = "data.frame", row.names = rownames, 
            dimensions = dimension.set))
    }, distances = {
        dclass <- matrix(0, nrow(newdata), J)
        for (i in seq(J)) dclass[, i] <- dist(newdata, means[i, ], 
            dimension)
        dimnames(dclass) <- list(dimnames(newdata)[[1]], dimnames(means)[[1]])
        return(dclass)
    })
}

