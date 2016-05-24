dynrq <- function (formula, tau = 0.5, data, subset, weights, na.action, method = "br", 
    contrasts = NULL, start = NULL, end = NULL, ...) 
{
    stopifnot(requireNamespace("zoo"))
    Zenv <- new.env(parent = environment(formula))
    assign("dynformula", function(x) structure(x, class = unique(c("dynformula", 
        oldClass(x)))), envir = Zenv)
    assign("L", function(x, k = 1) {
        if (length(k) > 1) {
            rval <- lapply(k, function(i) lag(x, k = -i))
            rval <- if (inherits(x, "ts")) 
                do.call("ts.intersect", rval)
            else do.call("zoo::merge", c(rval, list(all = FALSE)))
            colnames(rval) <- k
        }
        else {
            rval <- lag(x, k = -k)
        }
        return(rval)
    }, envir = Zenv)
    assign("d", function(x, lag = 1) diff(x, lag = lag), envir = Zenv)
    assign("season", function(x, ref = NULL) {
        freq <- frequency(x)
        stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), 
            TRUE))
        freq <- ofreq <- round(freq)
        freq <- if (freq == 12) 
            month.abb
        else if (freq == 4) 
            paste("Q", 1:4, sep = "")
        else 1:freq
        rval <- factor(zoo::coredata(cycle(x)), labels = freq)
        if (!is.null(ref)) 
            rval <- relevel(rval, ref = ref)
        rval <- zoo::zoo(rval, zoo::index(x), ofreq)
        return(rval)
    }, envir = Zenv)
    assign("trend", function(x, scale = TRUE) {
        freq <- ofreq <- if (inherits(x, "ts")) 
            frequency(x)
        else attr(x, "frequency")
        if (is.null(freq) | !scale) 
            freq <- 1
        stopifnot(freq >= 1 && identical(all.equal(freq, round(freq)), 
            TRUE))
        freq <- round(freq)
        rval <- zoo::zoo(seq_along(zoo::index(x))/freq, zoo::index(x), frequency = ofreq)
        return(rval)
    }, envir = Zenv)
    assign("harmon", function(x, order = 1) {
        freq <- frequency(x)
        stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), 
            TRUE))
        freq <- round(freq)
        order <- round(order)
        stopifnot(order <= freq/2)
        rval <- outer(2 * pi * zoo::index(x), 1:order)
        rval <- cbind(apply(rval, 2, cos), apply(rval, 2, sin))
        colnames(rval) <- if (order == 1) {
            c("cos", "sin")
        }
        else {
            c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, 
                sep = ""))
        }
        if ((2 * order) == freq) 
            rval <- rval[, -(2 * order)]
        return(rval)
    }, envir = Zenv)
    assign("model.frame.dynformula", function(formula, data = NULL, 
        subset = NULL, na.action = na.omit, drop.unused.levels = FALSE, 
        xlev = NULL, ...) {
        if (is.null(data)) 
            data <- parent.frame()
        if (!is.list(data)) 
            data <- as.list(data)
        args <- as.list(attr(terms(formula), "variables"))[-1]
        args$retclass <- "list"
        args$all <- FALSE
        formula <- terms(formula)
        attr(formula, "predvars") <- as.call(append(zoo::merge.zoo, args))
        attr(formula, "predvars")[[1]] <- as.name("merge.zoo")
        NextMethod("model.frame", formula = formula)
    }, envir = Zenv)
    if (missing(data)) 
        data <- Zenv
    orig.class <- if (is.data.frame(data) || is.environment(data)) 
        class(eval(attr(terms(formula), "variables")[[2]], data, 
            Zenv))
    else class(data)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- as.call(list(as.name("dynformula"), mf[[2]]))
    mf <- eval(mf, envir = Zenv)
    mfna <- attr(mf, "na.action")
    if (length(zoo::index(mf[, 1])) > nrow(mf)) {
        for (i in 1:NCOL(mf)) attr(mf[, i], "index") <- attr(mf[, 
            i], "index")[-as.vector(mfna)]
    }
    is.zoofactor <- function(x) !is.null(attr(x, "oclass")) && 
        attr(x, "oclass") == "factor"
    for (i in 1:NCOL(mf)) if (is.zoofactor(mf[, i])) 
        mf[, i] <- zoo::coredata(mf[, i])
    mf1 <- mf[, 1]
    start <- if (is.null(start)) 
        1
    else {
        if (length(start) > 1) 
            start <- start[1] + (start[2] - 1)/frequency(mf1)
        start <- min(which(zoo::index(mf1) >= start))
    }
    end <- if (is.null(end)) 
        length(mf1)
    else {
        if (length(end) > 1) 
            end <- end[1] + (end[2] - 1)/frequency(mf1)
        end <- max(which(zoo::index(mf1) <= end))
    }
    if (end < start) {
        warning("empty model frame specified")
        mf <- head(mf, 0)
        mf1 <- head(mf1, 0)
    }
    else {
        mf <- mf[start:end, , drop = FALSE]
        mf1 <- mf1[start:end]
        if (!is.null(mfna)) 
            attr(mf, "na.action") <- mfna[as.vector(mfna) >= 
                start & as.vector(mfna) <= end]
    }
    if ("ts" %in% orig.class && zoo::is.regular(mf1, strict = TRUE)) {
        for (i in 1:ncol(mf)) if (!is.factor(mf[, i])) 
            mf[, i] <- as.ts(mf[, i])
    }
    if (all(orig.class == "numeric")) {
        for (i in 1:ncol(mf)) if (!is.factor(mf[, i])) 
            mf[, i] <- as.vector(mf[, i])
    }
    rownames(mf) <- zoo::index2char(zoo::index(mf1), frequency(mf1))
    mt <- attr(mf, "terms")
    attr(mt, "predvars") <- NULL
    attr(mt, "dataClasses") <- NULL
    Y <- model.response(mf, "numeric")
    weights <- as.vector(model.weights(mf))
    if (is.empty.model(mt)) {
        X <- NULL
        rval <- list(coefficients = numeric(0), residuals = Y, 
            fitted.values = 0 * Y, weights = weights, rank = 0, df.residual = length(Y))
    }
    else {
        Rho <- function(u, tau) u * (tau - (u < 0))
        eps <- .Machine$double.eps^(2/3)
        X <- model.matrix(mt, mf, contrasts)
        if (length(tau) > 1) {
            if (any(tau < -eps) || any(tau > 1 + eps))
                stop("invalid tau:  taus should be >= 0 and <= 1")
            coef <- matrix(0, ncol(X), length(tau))
            rho <- rep(0, length(tau))
            fitted <- resid <- matrix(0, nrow(X), length(tau))
            for (i in 1:length(tau)) {
                z <- {
                    if (length(weights))
                      rq.wfit(X, Y, tau = tau[i], weights, method, ...)
                    else rq.fit(X, Y, tau = tau[i], method, ...)
                }
                coef[, i] <- z$coefficients
                resid[, i] <- z$residuals
                rho[i] <- sum(Rho(z$residuals, tau[i]))
                fitted[, i] <- Y - z$residuals
            }
            taulabs <- paste("tau=", format(round(tau, 3)))
            dimnames(coef) <- list(dimnames(X)[[2]], taulabs)
            dimnames(resid) <- list(dimnames(X)[[1]], taulabs)
            rval <- z
            rval$coefficients <- coef
            rval$residuals <- resid
            rval$fitted.values <- fitted
            class(rval) <- c("dynrqs", "rqs")
        }
    else {
        rval <- {
            if(length(weights))
                rq.wfit(X, Y, tau = tau, weights, method, ...)
            else rq.fit(X, Y, tau = tau, method, ...)
        }
        dimnames(rval$residuals) <- list(dimnames(X)[[1]], NULL)
        rho <- sum(Rho(rval$residuals, tau))
        class(rval) <- "rq"
        class(rval) <- c("dynrq", "rq")
        }

    }
    rval$na.action <- attr(mf, "na.action")
    rval$contrasts <- attr(X, "contrasts")
    rval$xlevels <- .getXlevels(mt, mf)
    rval$call <- cl
    rval$tau <- tau
    rval$terms <- mt
    rval$model <- mf
    rval$index <- zoo::index(mf1)
    rval$frequency <- frequency(mf1)
    rval$residuals <- drop(rval$residuals)
    rval$X <- X
    rval$y <- Y
    rval$rho <- rho
    rval$method <- method
    rval$fitted.values <- drop(rval$fitted.values)

    return(rval)
}

index.dynrq <- function(x, ...) {
  x$index
}

start.dynrq <- function(x, ...) {
  start(x$residuals)
}

end.dynrq <- function(x, ...) {
  end(x$residuals)
}

print.dynrq <- function(x, ...) {
  rx <- residuals(x)
  cat(paste("\nDynamic quantile regression \"", class(rx)[1], "\" data:\n", sep = ""))
  cat(paste("Start = ", zoo::index2char(zoo::index(rx)[1], x$frequency),
            ", End = ",   zoo::index2char(zoo::index(rx)[length(rx)], x$frequency), "\n", sep = ""))
  NextMethod()
}
print.dynrqs <- function(x, ...) {
  rx <- residuals(x)
  ix <- dimnames(rx)[[1]]
  cat(paste("\nDynamic quantile regression \"", class(rx)[1], "\" data:\n", sep = ""))
  cat(paste("Start = ", zoo::index2char(ix[1], x$frequency),
            ", End = ",   zoo::index2char(ix[length(ix)], x$frequency), "\n", sep = ""))
  NextMethod()
}

summary.dynrqs <- function(object, vcov. = NULL, df = NULL, ...) {
  rval <- NextMethod()
  #rval$frequency <- object$frequency
  #class(rval) <- c("summary.dynrqs", class(rval))
  return(rval)
}

print.summary.dynrq <- function(x, ...) {
  rx <- residuals(x)
  x$residuals <- zoo::coredata(x$residuals)
  cat(paste("\nDynamic quantile  regression \"", class(rx)[1], "\" data:\n", sep = ""))
  cat(paste("Start = ", zoo::index2char(zoo::index(rx)[1], x$frequency),
            ", End = ",   zoo::index2char(zoo::index(rx)[length(rx)], x$frequency), "\n", sep = ""))
  NextMethod()
}
print.summary.dynrqs <- function(x, ...) {
  lapply(x, print)
}
