## compute the squared residual length statistic for a constrained
## ordination given a single constraint

`residLen` <- function(X, env, passive,
                       method = c("cca", "rda")) {
    ## merge X and passive
    dat <- join(X, passive, type = "left") ## Think this should have type = "left" ?
    X <- dat[[1]]
    passive <- dat[[2]]
    ## check env is same length as nrow(X)
    if(!isTRUE(all.equal(length(env), nrow(X))))
        stop("'X' and 'env' imply different numbers of observations")
    ## ordinate
    if(missing(method))
        method <- "cca"
    method <- match.arg(method)
    FUN <- match.fun(method)
    ## ordinate
    ord <- FUN(X = X, Y = env)
    ## fitted values
    fit <- fitted(ord, type = "response")
    ## colSums of X
    ypk <- colSums(X)
    ## predict locations for the passive
    pred <- fittedY(ord, passive, ypk)
    ## calc sqdists
    if(method == "cca") {
        train <- sqrlUnimodal(X, ypk, fit)
        passive <- sqrlUnimodal(passive, ypk, pred)
    } else {
        train <- sqrlLinear(X, fit)
        passive <- sqrlLinear(passive, pred)
    }
    ## system.call
    .call <- match.call()
    .call[[1]] <- as.name("residLen")
    ## residual lengths for X
    res <- list(train = train, passive = passive,
                ordination = ord, call = .call)
    class(res) <- "residLen"
    attr(res, "method") <- method
    return(res)
}

fittedY <- function(ord, newdata, colsum) {
    ## Fitted values of response for samples
    ## Arguments:
    ## ord     = vegan ordination object (CCA/RDA)
    ## newdata = matrix of passive species data with same
    ##           columns as that used to fit ord
    ## colsum  = colum sums for the training
    ## species scores
    b <- predict(ord, type = "sp")
    ## site scores
    xi <- predict(ord, newdata = newdata, type = "wa")
    ## predict fitted values
    if(inherits(ord, "rda")) {
        fik <- xi %*% t(b)
    } else {
        fik <- sweep(1 + (xi %*% t(b)), 2,
                     colsum / sum(colsum), "*")
        fik <- sweep(fik, 1, rowSums(newdata), "*")
    }
    ## fitted values
    return(fik)
}

sqrlUnimodal <- function(Y, colsum, fitted) {
    ## Squared residual length for unimodal methods
    ## Arguments:
    ## Y      = species data (training or passive) to which
    ##          we want to find the residual length
    ## colsum = column sums for the training species data
    ## fitted = fitted species data for either the training
    ##          or passive samples
    yip <- rowSums(Y)
    ypk <- colsum
    ypp <- sum(ypk)
    A <- (sweep(Y, 2, ypk, "/") - sweep(fitted, 2, ypk, "/"))^2
    B <- sweep(A, 2, ypk, "*")
    C <- rowSums(B / ypp)
    res <- (1 / (yip / ypp))^2 * C
    return(res)
}

sqrlLinear <- function(Y, fitted) {
    ## Squared residual length for linear methods
    ## Arguments:
    ## Y      = species data (training or passive) to which
    ##          we want to find the residual length
    ## fitted = fitted species data for either the training
    ##          or passive samples
    res <- rowSums((Y - fitted)^2) / ncol(Y)
    return(res)
}
