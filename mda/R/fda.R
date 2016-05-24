fda <-
function (formula = formula(data), data = sys.frame(sys.parent()), 
    weights, theta, dimension = J - 1, eps = .Machine$double.eps, 
    method = polyreg, keep.fitted = (n * dimension < 5000), ...) 
{
    this.call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m[[1]] <- as.name("model.frame")
    m <- m[match(names(m), c("", "formula", "data", "weights"), 
        0)]
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    g <- model.extract(m, "response")
#    attr(Terms, "intercept") <- 0 This fails if a factor is in the model formula
    x <- model.matrix(Terms, m)
    if(attr(Terms, "intercept"))x=x[,-1,drop=FALSE]
    dd <- dim(x)
    n <- dd[1]
    weights <- model.extract(m, weights)
    if (!length(weights)) 
        weights <- rep(1, n)
    else if (any(weights < 0)) 
        stop("negative weights not allowed")
    if (length(g) != n) 
        stop("g should have length nrow(x)")
    fg <- factor(g)
    prior <- table(fg)
    prior <- prior/sum(prior)
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    iswt <- FALSE
    if (missing(weights)) 
        dp <- table(g)/n
    else {
        weights <- (n * weights)/sum(weights)
        dp <- tapply(weights, g, sum)/n
        iswt <- TRUE
    }
    if (missing(theta)) 
        theta <- contr.helmert(J)
    theta <- contr.fda(dp, theta)
    Theta <- theta[g, , drop = FALSE]
    fit <- method(x, Theta, weights, ...)
    if (iswt) 
        Theta <- Theta * weights
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = fit, call = this.call), 
                         class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(cnames, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    obj <- structure(list(percent.explained = pe, values = lambda, 
        means = means, theta.mod = thetan, dimension = dimension, 
        prior = prior, fit = fit, call = this.call, terms = Terms), 
        class = "fda")
    obj$confusion <- confusion(predict(obj), fg)
    if (!keep.fitted) 
        obj$fit$fitted.values <- NULL
    obj
}

