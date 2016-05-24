loclda <- function(x, ...)
    UseMethod("loclda")

loclda.default <- function(x, grouping, weight.func = function(x) 1/exp(x),
                           k = nrow(x), weighted.apriori = TRUE, ...)
{
    cl <- match.call()
    cl[[1]] <- as.name("loclda")
    structure(list(learn = x, grouping = grouping, lev = levels(grouping),
                   weight.func = weight.func, k = k, weighted.apriori = weighted.apriori,
                   call = cl),
        class = "loclda")
}

loclda.formula <- function(formula, data = NULL, ...,
                           subset, na.action = na.fail)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0)
        x <- x[, -xint, drop = FALSE]
    res <- loclda(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("loclda")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- .getXlevels(Terms, m)
    res$na.action <- attr(m, "na.action")
    res
}

loclda.matrix <- function (x, grouping, ..., subset, na.action = na.fail)
{
    if (!missing(subset)){
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if (!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x),
            class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- loclda.default(x, grouping, ...)
    cl <- match.call()
    cl[[1]] <- as.name("loclda")
    res$call <- cl
    res
}



loclda.data.frame <- function (x, ...)
{
    res <- loclda.matrix(structure(data.matrix(x), class = "matrix"),
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("loclda")
    res$call <- cl
    res
}


predict.loclda <- function(object, newdata,...)
{
    weighted <- function (z, Z, y, weight.func, k, weighted.apriori)
    {
        unique.y <- sort(unique(y))
        Zn <- kronecker(rep(1, nrow(Z)), t(z)) - Z
        distances <- sqrt(rowSums(Zn^2))
        Knn <- sort(distances)[k]
        indic <- distances <= Knn
        kdist <- distances[indic] / Knn
        weighted <- sapply(kdist, weight.func)
        y <- y[indic]
        Z <- Z[indic, , drop = FALSE]
        indic2 <- y %in% (unique.y[table(y) > 1])
        y <- y[indic2]
        weighted <- weighted[indic2]
        if (is.vector(Z))
            Z <- matrix(Z, ncol = 1)
        Z <- Z[indic2, , drop = FALSE]
        apriori.table <- sapply(unique.y, function(x) sum(weighted[y == x]))
        apriori.table <- apriori.table / sum(weighted)
        if (!weighted.apriori)
        {
            indic3 <- apriori.table != 0
            apriori.table [indic3] <- 1/sum (indic3)
        }
        return(list(y = y, Weights = weighted, X = Z, apriori = apriori.table))
    }

    cov.pooled <- function (X, y, Weights)
    {
        y <- as.numeric(y)
        hilfe <- unique(y)
        if (is.vector(X))
            X <- matrix(X, ncol = 1)
        cov.vecs <- sapply(hilfe, function(r){
            temp <- y == r
            cov.wt(matrix(X[temp, ], ncol = ncol(X)), wt = Weights[temp])$cov
        })
        if(ncol(X) == 1)
        cov.vecs <- matrix(cov.vecs, nrow = 1)
        ni <- sapply(hilfe, function(r) sum(y == r))
        cov.vecs <- cov.vecs %*% diag(ni, nrow = length(ni))
        1 / (nrow(X) - length(hilfe)) * matrix(rowSums(cov.vecs), nrow = ncol(X))
    }

    mu.weighted <- function (X, unique.y, y2, Weights, i)
    {
        if (is.vector(X))
            X <- matrix(X, ncol = 1)
        temp <- y2 == unique.y[i]
        Xnow <- X[temp, , drop = FALSE]
        Wnow <- Weights[temp]
        as.vector(t(Wnow) %*% Xnow / sum(Wnow))
    }

    ML <- function (x, X, y, y2, Weights, apriori)
    {
        unique.y <- sort(unique(y))
        sigma <- cov.pooled(X = X, y = y2, Weights = Weights)
        sigma <- solve(sigma, tol = 1e-100)
        mu.distances <- rep (Inf, length (unique.y))
        posterior.loglik <- -mu.distances
        for (i in seq(along=posterior.loglik)[apriori!= 0])
        {
            mu <- mu.weighted(X = X, unique.y = unique.y, y2 = y2, Weights = Weights, i = i)
            temp <- x - mu
            mu.distances[i] <- drop (crossprod (temp))
            posterior.loglik[i] <-(-1/2) * t(temp) %*% sigma %*% (temp) + log (apriori[i])
        }
        names (posterior.loglik) <- names (mu.distances) <- as.character(unique.y)
        c (mu.distances, posterior.loglik)
    }

    llda <- function (x, learn, y, weight.func, k, weighted.apriori)
    {
        Weighted <- weighted(z = x, Z = learn, y = y,
                             weight.func = weight.func, k = k, weighted.apriori = weighted.apriori)
        ML(x = x, X = Weighted$X, y = y, y2 = Weighted$y, Weights = Weighted$Weights,
           apriori = Weighted$apriori)
    }

    if (!inherits(object, "loclda"))
            stop("object not of class", " 'loclda'")
        if (!is.null(Terms <- object$terms)) {
            if (missing(newdata))
                newdata <- model.frame(object)
            else {
                newdata <- model.frame(as.formula(delete.response(Terms)),
                    newdata, na.action = function(x) x, xlev = object$xlevels)
            }
            x <- model.matrix(delete.response(Terms), newdata,
                              contrasts = object$contrasts)
            xint <- match("(Intercept)", colnames(x), nomatch = 0)
            if (xint > 0)
                x <- x[, -xint, drop = FALSE]
        }
        else {
            if (missing(newdata)){
                if (!is.null(sub <- object$call$subset))
                    newdataa <- eval.parent(parse(
                        text = paste(deparse(object$call$x, backtick = TRUE),
                                     "[", deparse(sub, backtick = TRUE), ",]")))
                else newdata <- eval.parent(object$call$x)
                if (!is.null(nas <- object$call$na.action))
                    newdata <- eval(call(nas, newdata))
            }
            if (is.null(dim(newdata)))
                dim(newdata) <- c(1, length(newdata))
            x <- as.matrix(newdata)
        }

    werte <- t(apply(x, 1, llda, learn = object$learn,
        y = object$grouping, weight.func = object$weight.func,
        k = object$k, weighted.apriori = object$weighted.apriori))
    object.levels <- object$lev
    indic4 <- seq (along = object.levels)
    mu.distances <- werte[,indic4]
    posterior.loglik <- werte[,-indic4]
    if (nrow (newdata) == 1) {
        mu.distances <- matrix (mu.distances, nrow = 1)
        posterior.loglik <- matrix (posterior.loglik, nrow = 1)
    }    
    posterior <- exp (posterior.loglik)
    classes <- factor(max.col(posterior), levels = indic4,
                      labels = object.levels)
    dim.posterior <- dim (posterior)
    all.zero <-  logical (dim.posterior[1])
    for (i in seq (along = classes)) 
        if (all (posterior.loglik[i,] < log(1e-150))) {
            classes[i] <- object.levels[which.min(mu.distances[i,])]
            all.zero[i] <- TRUE
        }
    all.zero <- which (all.zero)
    posterior <- apply (posterior, 1, function (x) 
        if(sum(x) != 0) 
            x <- x/sum(x) 
        else 
            x <- rep (1/dim.posterior[2], dim.posterior[2])
    )
    result <- list(class = classes, posterior = t(posterior), all.zero = all.zero)
    return(result)
}

print.loclda <- function (x, ...)
{
    cat("Call:\n")
    print(x$call, ...)
    cat("\nWeighting function:\n")
    print (x$weight.func, ...)
    cat("\nNumber of next neighbours that will be used for prediction:\n")
    print(x$k, ...)
    cat("\nUsage of weighted a priori probabilities:\n")
    print(x$weighted.apriori, ...)
    invisible(x)
}
