mda.start <-
function (x, g, subclasses = 3, trace.mda.start = FALSE, start.method = c("kmeans", 
    "lvq"), tries = 5, criterion = c("misclassification", "deviance"), 
    ...) 
{
    start.method <- match.arg(start.method)
    criterion <- match.arg(criterion)
    name.criterion <- switch(criterion, misclassification = "Misclassification Error", 
        deviance = "Deviance(multinomial)")
    starter <- get(paste(start.method, "start", sep = "."), mode = "function")
    fg <- factor(g)
    cnames <- levels(fg)
    prior <- table(fg)
    J <- length(cnames)
    n <- length(g)
    g <- as.numeric(fg)
    best.ll <- 1/.Machine$double.eps
    for (try in seq(tries)) {
        start <- starter(x, fg, subclasses)
        weights <- start$weights
        if (criterion == "misclassification") {
            pg <- lvqtest(start, x)
            ll <- attr(confusion(pg, fg), "error")
        }
        else {
            subclasses <- sapply(weights, ncol)
            Assign <- split(seq(sum(subclasses)), rep(seq(J), 
                subclasses))
            names(Assign) <- cnames
            fit <- mda.fit(x, g, weights, assign.theta = Assign, 
                Rj = subclasses, eps = .Machine$double.eps, method = polyreg, 
                ...)
            dmat <- exp(-0.5 * predict.fda(fit, type = "distance"))
            sub.prior <- fit$sub.prior
            pclass <- matrix(1, n, J)
            for (j in seq(J)) {
                priorj <- sub.prior[[j]]
                ass <- Assign[[j]]
                TT <- dmat[, ass, drop = FALSE]
                TT <- TT * outer(rep(1, n), priorj)
                TTot <- drop(TT %*% rep(1, length(ass)))
                wmj <- TT[g == j, , drop = FALSE]/TTot[g == j]
                pclass[, j] <- prior[j] * TTot
                dimnames(wmj) <- list(NULL, paste("s", seq(along = ass), 
                  sep = ""))
                weights[[j]] <- wmj
            }
            pclass <- pclass/drop(pclass %*% rep(1, J))
            ll <- llmult(pclass, g)
        }
        if (trace.mda.start) 
            cat(start.method, "start   \t", name.criterion, format(round(ll, 
                5)), "\n")
        if (ll < best.ll) {
            keep.weights <- weights
            best.ll <- ll
        }
    }
    structure(keep.weights, criterion = best.ll, name.criterion = name.criterion)
}

