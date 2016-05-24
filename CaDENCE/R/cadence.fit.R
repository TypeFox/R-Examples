cadence.fit <-
function (x, y, iter.max = 500, n.hidden = 2, hidden.fcn = tanh, 
          distribution = NULL, sd.norm = Inf, init.range = c(-0.5, 0.5),
          method = c("optim", "psoptim", "Rprop"), n.trials = 1,
          trace = 0, maxit.Nelder = 2000, trace.Nelder = 0,
          swarm.size = NULL, vectorize = TRUE,
          delta.0 = 0.1, delta.min = 1e-06, delta.max = 50, epsilon = 1e-08,
          range.mult = 2, step.tol = 1e-08, f.target = -Inf,
          f.cost = cadence.cost, max.exceptions = 500)
{
    if (!is.matrix(x)) 
        stop("\"x\" must be a matrix")
    if (!is.matrix(y)) 
        stop("\"y\" must be a matrix")
    if (any(c(n.hidden, n.trials) <= 0)) 
        stop("invalid \"n.hidden\", or \"n.trials\"")
    if (is.null(distribution)) {
        distribution <- list(density.fcn = dnorm, parameters = c("mean", 
            "sd"), output.fcns = c(identity, exp))
        warning("unspecified distribution; defaulting to dnorm")
    }
    if (identical(hidden.fcn, identity)) 
        n.hidden <- length(distribution$parameters)
    stationary <- FALSE
    if (sum(distribution$parameters %in% distribution$parameters.fixed) == 
        length(distribution$parameters)) {
        x <- x[, 1, drop = FALSE] * 0
        hidden.fcn <- identity
        n.hidden <- 1
        stationary <- TRUE
        warning("fitting stationary model")
    }
    x <- scale(x)
    attr(x, "scaled:scale")[attr(x, "scaled:scale") == 0] <- 1
    x[is.nan(x)] <- 0
    method <- match.arg(method)
    fprime <- function(w.init, f, epsilon, f.init, ...) {
        if (is.null(f.init)) 
            f.init <- f(w.init, ...)
        gradient <- w.init * 0
        for (i in 1:length(w.init)) {
            w.plus <- w.init
            w.plus[i] <- w.plus[i] + epsilon
            f.plus <- f(w.plus, ...)
            gradient[i] <- (f.plus - f.init)/epsilon
        }
        gradient
    }
    fit <- list()
    for (nh in n.hidden) {
        NLL <- Inf
        cat("n.hidden =", nh, "--> ")
        for (i in 1:n.trials) {
            cat(i, "")
            exception <- TRUE
            n.exceptions <- 0
            while (exception) {
                weights <- cadence.initialize(x, nh, init.range, 
                  distribution)
                gradient <- fprime(weights, f.cost, epsilon, 
                  NULL, x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn, 
                  distribution = distribution, sd.norm = sd.norm, 
                  valid = rep(TRUE, length(weights)))
                valid.cur <- gradient^2 > epsilon
                weights <- weights[valid.cur]
                if (trace > 0) 
                  cat("nonzero weights =", valid.cur, "\n")
                if (method == "optim") {
                  output.cdn.cur <- try(suppressWarnings(optim(weights, 
                    f.cost, method = "Nelder",
                    control = list(maxit = maxit.Nelder, trace = trace.Nelder),
                    x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn,
                    distribution = distribution, sd.norm = sd.norm,
                    valid = valid.cur)), silent = trace == 0)
                  weights <- try(output.cdn.cur$par, silent = trace == 0)
                  if(iter.max > 0){
                      output.cdn.cur <- try(suppressWarnings(optim(weights, 
                        f.cost, method = "BFGS",
                        control = list(maxit = iter.max, reltol = step.tol,
                        trace = trace, REPORT = ifelse(trace<=0, 1, trace)),
                        x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn,
                        distribution = distribution, sd.norm = sd.norm,
                        valid = valid.cur)), silent = trace == 0)
                   }
                }
                else if (method == "psoptim") {
                  output.cdn.cur <- try(suppressWarnings(optim(weights, 
                    f.cost, method = "Nelder",
                    control = list(maxit = maxit.Nelder, trace = trace.Nelder),
                    x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn,
                    distribution = distribution, sd.norm = sd.norm,
                    valid = valid.cur)), silent = trace == 0)
                  weights <- try(output.cdn.cur$par, silent = trace == 0)
                  if(iter.max > 0){
                    w.lower <- weights - range.mult*diff(range(weights))
                    w.upper <- weights + range.mult*diff(range(weights))
                    if(is.null(swarm.size))
                        swarm.size <- floor(10+2*sqrt(length(weights)))
                    output.cdn.cur <- try(suppressWarnings(pso::psoptim(weights, 
                        f.cost, lower = w.lower, upper = w.upper,
                        control = list(maxit = iter.max, abstol = f.target,
                        vectorize = vectorize, s = swarm.size,
                        trace = trace, REPORT = ifelse(trace<=0, 1, trace)),
                        x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn,
                        distribution = distribution, sd.norm = sd.norm,
                        valid = valid.cur)), silent = trace == 0)
                   }
                }
                else if (method == "Rprop") {
                  output.cdn.cur <- try(suppressWarnings(optim(weights, 
                    f.cost, method = "Nelder",
                    control = list(maxit = maxit.Nelder, trace = trace.Nelder),
                    x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn,
                    distribution = distribution, sd.norm = sd.norm,
                    valid = valid.cur)), silent = trace == 0)
                  weights <- try(output.cdn.cur$par, silent = trace == 0)
                  if(iter.max > 0){
                    output.cdn.cur <- try(rprop(w = weights, f = f.cost, 
                        iterlim = iter.max, print.level = trace, 
                        delta.0 = delta.0, delta.min = delta.min, 
                        delta.max = delta.max, epsilon = epsilon, 
                        step.tol = step.tol, f.target = f.target, 
                        x = x, y = y, n.hidden = nh, hidden.fcn = hidden.fcn, 
                        distribution = distribution, sd.norm = sd.norm, 
                        valid = valid.cur), silent = trace == 0)
                    }
                }
                if (!class(output.cdn.cur) == "try-error") {
                  exception <- FALSE
                }
                else {
                  n.exceptions <- n.exceptions + 1
                  if (n.exceptions > max.exceptions) 
                    stop("max. number of exceptions reached")
                }
            }
            NLL.cur <- output.cdn.cur$value
            if (NLL.cur < NLL) {
                NLL <- NLL.cur
                output.cdn <- output.cdn.cur
                valid <- valid.cur
            }
        }
        weights <- output.cdn$par
        NLL <- f.cost(weights, x, y, nh, hidden.fcn, distribution, 
            sd.norm, valid)
        penalty <- attr(NLL, "penalty")
        attr(NLL, "penalty") <- NULL
        NLL <- NLL - penalty
        cat("* NLL =", NLL, "; penalty =", penalty)
        n.parms <- length(distribution$parameters)
        if (stationary) {
            k <- n.parms
        }
        else {
            if (identical(hidden.fcn, identity)) {
                k <- (n.parms - length(distribution$parameters.fixed)) * 
                  (ncol(x) + 1) + length(distribution$parameters.fixed)
            }
            else {
                k <- length(weights)
                if (!is.null(distribution$parameters.fixed)) 
                  k <- k - nh * length(distribution$parameters.fixed)
            }
        }
        n <- nrow(y)
        BIC <- 2 * NLL + k * log(n)
        AIC <- 2 * NLL + 2 * k
        AICc <- AIC + (2 * k * (k + 1))/(n - k - 1)
        weights.final <- valid * 0
        weights.final[valid] <- weights
        w <- cadence.reshape(x, weights.final, nh, distribution)
        attr(w, "n.hidden") <- nh
        attr(w, "k") <- k
        attr(w, "NLL") <- NLL
        attr(w, "penalty") <- penalty
        attr(w, "BIC") <- BIC
        attr(w, "AICc") <- AICc
        attr(w, "AIC") <- AIC
        if (sd.norm == Inf) {
            cat("; BIC =", BIC, "; AICc =", AICc, "; AIC =", 
                AIC)
        }
        cat("\n")
        attr(w, "hidden.fcn") <- hidden.fcn
        attr(w, "distribution") <- distribution
        attr(w, "x.center") <- attr(x, "scaled:center")
        attr(w, "x.scale") <- attr(x, "scaled:scale")
        attr(w, "stationary") <- stationary
        fit[[as.character(nh)]] <- w
    }
    if (sd.norm < Inf) {
        warning("The value of \"sd.norm\" < Inf. Treat values of \"k\", \"BIC\", \"AICc\", and \"AIC\" with caution.")
    }
    fit
}
