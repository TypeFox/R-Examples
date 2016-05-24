heavyPS <-
function(x, y, family = Student(df = 4), nseg = 20, deg = 3, ord = 2, lambda = 1, method = c("GCV", "none"), ngrid = 200, control = heavy.control())
{
    ## local functions, from Eilers and Marx (1996, 2010)
    tpower <- function(x, s, p) {
        # truncated p-th power function
        (x - s) ^ p * (x > s)
    }
    
    bbase <- function(x, xmin = min(x), xmax = max(x), nseg = 20, deg = 3) {
        # construct a B-spline basis of degree 'deg'
        dx <- (xmax - xmin) / nseg
        knots <- seq(xmin - deg * dx, xmax + deg * dx, by = dx)
        P <- outer(x, knots, tpower, deg)
        n <- dim(P)[2]
        K <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
        B <- (-1) ^ (deg + 1) * P %*% t(K)
        list(B = B, knots = knots)
    }
    
    penalty <- function(p, ord) {
        # construct a penalty matrix of order 'ord'
        K <- diff(diag(p), differences = ord)
        K
    }
    
    ## validating arguments
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    if (!is.numeric(y)) stop("'y' must be a numeric vector")
    ## remove all NAs
    OK <- complete.cases(x, y)
    x <- x[OK]
    y <- y[OK]
    
    ## initial computations
    cl <- match.call()
    xmin <- min(x)
    xmax <- max(x)
    z <- bbase(x, xmin = xmin, xmax = xmax, nseg = nseg, deg = deg)
    B <- z$B
    n <- nrow(B)
    p <- ncol(B)
    K <- penalty(p, ord)
    des <- list(B = B, K = K, knots = z$knots)
    dims <- c(n, p, deg, ord)
  
    ## extract family info
    if (!inherits(family, "heavy.family"))
        stop("Use only with 'heavy.family' objects")
    if (is.null(family$family))
        stop("'family' not recognized")
    kind <- family$which
    if ((kind < 0) || (kind > 4))
        stop("not valid 'family' object")
    settings <- c(kind, family$npars, unlist(family$pars))

    ## set control values and method
    if (missing(control))
        control <- heavy.control()
    ctrl <- unlist(control)[c(1:4, 6)]
    ctrl <- c(ctrl, 0)
    method <- match.arg(method)
  
    ## initial estimates
    aug <- rep(0, p - ord)
    z <- lsfit(rbind(B, sqrt(lambda) * K), c(y, aug), intercept = FALSE)
    a <- z$coef
    resid <- z$resid[1:n]
    pen <- sum(c(K %*% a)^2)
    scale <- (sum(resid^2) + lambda * pen) / n
    distances <- resid^2 / scale
    weights <- rep(1, n)
  
    ## set the storage mode
    storage.mode(y) <- "double"
    storage.mode(B) <- "double"
    storage.mode(K) <- "double"
    storage.mode(dims) <- "integer"

    ## call fitter subroutine
    now <- proc.time()
    z <- switch(method,
                "GCV"  = .C("ps_combined",
                            y = y, B = B, K = K, dims = dims,
                            settings = as.double(settings),
                            coef = as.double(a),
                            scale = as.double(scale),
                            lambda = as.double(lambda),
                            edf = as.double(0),
                            gcv = as.double(0),
                            pen = as.double(pen),
                            fitted = double(n),
                            resid = as.double(resid),
                            distances = as.double(distances),
                            weights = as.double(weights),
                            plogLik = as.double(0),
                            control = as.double(ctrl)),
                "none" = .C("ps_fit",
                            y = y, B = B, K = K, dims = dims,
                            settings = as.double(settings),
                            coef = as.double(a),
                            scale = as.double(scale),
                            lambda = as.double(lambda),
                            edf = as.double(0),
                            gcv = as.double(0),
                            pen = as.double(pen),
                            fitted = double(n),
                            resid = as.double(resid),
                            distances = as.double(distances),
                            weights = as.double(weights),
                            plogLik = as.double(0),
                            control = as.double(ctrl)),
                stop(paste("unimplemented method:", method))
                )
    speed <- proc.time() - now

    ## compute curve on grid
    xgrid <- seq(xmin, xmax, length = ngrid)
    B <- bbase(xgrid, xmin, xmax, nseg = nseg, deg = deg)$B
    ygrid <- c(B %*% z$coef)
  
    ## creating the output object
    out <- list(call = cl,
                design = des,
                dims = dims,
                xmin = xmin,
                xmax = xmax,
                nseg = nseg,
                method = method,
                family = family,
                settings = z$settings,
                coefficients = z$coef,
                scale = z$scale,
                lambda = z$lambda,
                fitted.values = z$fitted,
                residuals = z$resid,
                plogLik = z$plogLik,
                edf = z$edf,
                gcv = exp(z$gcv),
                pen = z$pen,
                numIter = z$control[6],
                control = z$control[1:5],
                weights = z$weights,
                distances = z$distances,
                xgrid = xgrid,
                ygrid = ygrid,
                speed = speed,
                converged = FALSE)
    if (!control$fix.shape) {
        if ((kind > 1) && (kind < 4)) {
            df <- signif(out$settings[3], 6)
            out$family$call <- call(out$family$family, df = df)
        }
    }
    if (out$numIter < control$maxIter)
        out$converged <- TRUE
    if (!control$fix.shape)
        out$shape <- z$settings[-c(1:2)]
    class(out) <- "heavyPS"
    out
}

print.heavyPS <-
function(x, digits = max(4, getOption("digits") - 4), ...)
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        cl$family <- x$family$call
        dput(cl, control = NULL)
    }
    if (x$converged)
        cat("Converged in", x$numIter, "iterations\n")
    else
        cat("Maximum number of iterations exceeded\n")
    cat("\nNumber of observations:", x$dims[1], "\n")
    cat("Equivalent number of parameters:", format(round(x$edf, 2)), "\n")
    cat("Residual scale estimate:", format(signif(x$scale, digits)), "\n")
    invisible(x)
}

summary.heavyPS <- function(object, ...)
{
    class(object) <- "summary.heavyPS"
    object
}

print.summary.heavyPS <- function(x, digits = max(4, getOption("digits") - 4), ...)
{
    cat("Penalized spline under heavy-tailed distributions\n")
    dnames <- paste(x$call$x, x$call$y, sep = " and ")
    cat(" Data:", paste(dnames, ";", sep = ""))
    print(x$family)
    cat("\nNumber of observations:", x$dims[1], "\n")
    cat("Equivalent number of parameters:", format(round(x$edf, 2)), "\n")
    cat("Residual scale estimate:", format(signif(x$scale, digits)), "\n")
    if (x$method == "EM")
        cat("Smoothing parameter estimate:", paste(format(signif(x$lambda, digits)), ",", sep = ""))
    else if (x$method == "GCV")
        cat("Smoothing parameter estimate:", paste(format(signif(x$lambda, digits)), ",", sep = ""))
    else if (x$method == "none")
        cat("Smoothing parameter:", paste(format(signif(x$lambda, digits)), ",", sep = ""))
    cat(" WGCV:", format(signif(x$gcv, digits)), "\n")
    cat("\nControl settings:\n")
    cat("  segments: ", x$nseg, "\n")
    cat("  degree  : ", x$dims[3], "\n")
    cat("  order   : ", x$dims[4])
    cat("\n")
    invisible(x)
}
