monmlp.fit <-
function(x, y, hidden1, hidden2=0, iter.max=5000, n.trials=1, n.ensemble=1,
         bag=FALSE, cases.specified=NULL, iter.stopped=NULL, scale.y=TRUE,
         Th=tansig, To=linear, Th.prime=tansig.prime, To.prime=linear.prime,
         monotone=NULL, init.weights=c(-0.5, 0.5), max.exceptions=10,
         silent=FALSE, ...)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    if (hidden1 <= 0 | hidden2 < 0) stop("invalid no. of hidden nodes")
    if (any(c(iter.max, n.trials, n.ensemble) <= 0))
        stop("invalid \"iter.max\", \"n.trials\", or \"n.ensemble\"")
    if (!is.null(cases.specified) & length(cases.specified) != n.ensemble)
        stop("invalid \"cases.specified\"")
    x.raw <- x
    y.raw <- y
    x <- scale(x.raw)
    attr(x, "scaled:scale")[attr(x, "scaled:scale")==0] <- 1
    x[is.nan(x)] <- 0
    if (scale.y) y <- scale(y.raw)
    cases <- 1:nrow(x)
    oob <- NULL
    w.ens <- list()
    for (ens in 1:n.ensemble){
        if (!silent) cat("** Ensemble", ens, "\n")
        if (bag){
            if (is.null(cases.specified)){
                cases <- sample(nrow(x), replace=TRUE)
            } else{
                cases <- cases.specified[[ens]]
            }
            oob <- which(!(1:nrow(x) %in% cases))
            if (!silent) cat("** Bagging on\n")
        }
        if (bag & !is.null(iter.stopped)){
            if (!silent) cat("** Stopped training on\n")
            cost.best <- Inf
            iter <- 0
            while (iter < iter.max){
                if (iter==0){
                    weights <- init.weights
                    n.trials.stopped <- n.trials
                } else {
                    n.trials.stopped <- 1
                }
                iter <- iter + iter.stopped
                fit.ens <- monmlp.nlm(x=x[cases,,drop=FALSE],
                                      y=y[cases,,drop=FALSE],
                                      hidden1=hidden1, hidden2=hidden2,
                                      iter.max=iter.stopped,
                                      n.trials=n.trials.stopped,
                                      Th=Th, To=To, Th.prime=Th.prime,
                                      To.prime=To.prime, monotone=monotone,
                                      init.weights=weights,
                                      max.exceptions=max.exceptions,
                                      silent=silent, ...)
                weights <- fit.ens$weights
                code <- fit.ens$code
                w <- list(monmlp.reshape(x=x, y=y, weights=weights,
                                         hidden1=hidden1, hidden2=hidden2))
                attr(w, "Th") <- Th
                attr(w, "To") <- To
                attr(w, "monotone") <- monotone
                attr(w, "x.center") <- rep(0, ncol(x))
                attr(w, "x.scale") <- rep(1, ncol(x))
                attr(w, "y.center") <- rep(0, ncol(y))
                attr(w, "y.scale") <- rep(1, ncol(y))
                pred <- monmlp.predict(x=x[oob,,drop=FALSE], weights=w)
                cost <- mean((y[oob,,drop=FALSE]-pred)^2)
                if (!silent) cat("\t  --->", iter, cost, "\n")
                if (cost < cost.best){
                    weights.best <- weights
                    cost.best <- cost
                    iter.best <- iter
                }
                if (code <= 2) break
            }
            if (!silent) cat("**", iter.best, cost.best, "\n\n")
            w <- monmlp.reshape(x=x, y=y, weights=weights.best, hidden1=hidden1,
                                hidden2=hidden2)
            attr(w, "iter.best") <- iter.best
            attr(w, "cost.best") <- cost.best
        } else{
            iter.stopped <- NULL
            fit.ens <- monmlp.nlm(x=x[cases,,drop=FALSE],
                                  y=y[cases,,drop=FALSE],
                                  hidden1=hidden1, hidden2=hidden2,
                                  iter.max=iter.max, n.trials=n.trials, Th=Th,
                                  To=To, Th.prime=Th.prime, To.prime=To.prime,
                                  monotone=monotone, init.weights=init.weights,
                                  max.exceptions=max.exceptions,
                                  silent=silent, ...)
            weights <- fit.ens$weights
            cost <- fit.ens$cost
            if (!silent) cat("**", cost, "\n\n")
            w <- monmlp.reshape(x=x, y=y, weights=weights, hidden1=hidden1,
                                hidden2=hidden2)
        }
        attr(w, "oob") <- oob
        w.ens[[ens]] <- w
    }
    attr(w.ens, "x") <- x.raw
    attr(w.ens, "y") <- y.raw
    attr(w.ens, "Th") <- Th
    attr(w.ens, "To") <- To
    attr(w.ens, "Th.prime") <- Th.prime
    attr(w.ens, "To.prime") <- To.prime
    attr(w.ens, "monotone") <- monotone
    attr(w.ens, "bag") <- bag
    attr(w.ens, "iter.max") <- iter.max
    attr(w.ens, "iter.stopped") <- iter.stopped
    attr(w.ens, "x.center") <- attr(x, "scaled:center")
    attr(w.ens, "x.scale") <- attr(x, "scaled:scale")
    if (scale.y){
        attr(w.ens, "y.center") <- attr(y, "scaled:center")
        attr(w.ens, "y.scale") <- attr(y, "scaled:scale")
    } else{
        attr(w.ens, "y.center") <- rep(0, ncol(y))
        attr(w.ens, "y.scale") <- rep(1, ncol(y))
    }
    y.pred <- monmlp.predict(x=x.raw, weights=w.ens)
    attr(w.ens, "y.pred") <- y.pred
    w.ens
}

