### generic method for asymptotic null distributions
setGeneric("AsymptNullDistribution",
    function(object, ...) {
        standardGeneric("AsymptNullDistribution")
    }
)

### method for scalar test statistics
setMethod("AsymptNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object, ...) {
        p <- function(q) pnorm(q)
        q <- function(p) qnorm(p)
        pvalue <- function(q)
            switch(object@alternative,
                "less"      = p(q),
                "greater"   = 1 - p(q),
                "two.sided" = 2 * min(p(q), 1 - p(q)))

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,
            q = q,
            d = function(x) dnorm(x),
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Univariate Normal Distribution")
    }
)

### method for max-type test statistics
setMethod("AsymptNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        corr <- cov2cor(covariance(object))
        pq <- length(expectation(object))
        p <- function(q, ..., conf.int = FALSE) {
            p <- function(q, ..., conf.int = conf.int)
                switch(object@alternative,
                    "less"      = pmvn(lower = q, upper = Inf,
                                       mean = rep.int(0, pq),
                                       corr = corr, ..., conf.int = conf.int),
                    "greater"   = pmvn(lower = -Inf, upper = q,
                                       mean = rep.int(0, pq),
                                       corr = corr, ..., conf.int = conf.int),
                    "two.sided" = pmvn(lower = -abs(q), upper = abs(q),
                                       mean = rep.int(0, pq),
                                       corr = corr, ..., conf.int = conf.int))
            if (length(q) > 1)
                vapply(q, p, NA_real_, conf.int = conf.int)
            else
                p(q, conf.int = conf.int)
        }
        q <- function(p, ...) {
            q <- function(p, ...)
                qmvn(p, mean = rep.int(0, pq), corr = corr, ...)
            if (length(p) > 1)
                vapply(p, q, NA_real_)
            else
                q(p)
        }
        pvalue <- function(q) {
            pvalue <- 1 - p(q, conf.int = TRUE)
            ci <- 1 - attr(pvalue, "conf.int")[2L:1L]
            attr(ci, "conf.level") <- attr(attr(pvalue, "conf.int"), "conf.level")
            attr(pvalue, "conf.int") <- ci
            class(pvalue) <- "MCp"
            pvalue
        }

        new("AsymptNullDistribution",
            seed = seed,
            p = p,
            q = q,
            d = function(x) NA,
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Multivariate Normal Distribution",
            parameters = list(corr = corr))
    }
)

### method for quad-type test statistics
setMethod("AsymptNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        p <- function(q) pchisq(q, df = object@df)
        q <- function(p) qchisq(p, df = object@df)
        pvalue <- function(q) 1 - p(q)

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,
            q = q,
            d = function(d) dchisq(d, df = object@df),
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Chi-Squared Distribution",
            parameters = list(df = object@df))
    }
)


### generic method for exact null distributions
setGeneric("ExactNullDistribution",
    function(object, ...) {
        standardGeneric("ExactNullDistribution")
    }
)

### method for scalar test statistics
setMethod("ExactNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object,
        algorithm = c("auto", "shift", "split-up"), ...) {
            algorithm <- match.arg(algorithm)
            if (object@paired) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for paired samples")
                int <- is_integer(object@ytrans[, 1L], ...)
                if (int)
                    SR_shift_1sample(object, fact = attr(int, "fact"))
                else
                    stop("cannot compute exact distribution with real-valued scores")
            } else if (is_2sample(object)) {
                if (algorithm == "split-up")
                    vdW_split_up_2sample(object)
                else {
                    int <- is_integer(object@ytrans[, 1L], ...)
                    if (int)
                        SR_shift_2sample(object, fact = attr(int, "fact"))
                    else if (algorithm == "auto")
                        vdW_split_up_2sample(object)
                    else
                        stop("cannot compute exact distribution with real-valued scores")
                }
            } else
                stop(sQuote("object"), " is not a two-sample problem")
        }
)

### method for quad-type test statistics
setMethod("ExactNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object,
        algorithm = c("auto", "shift", "split-up"), ...) {
            algorithm <- match.arg(algorithm)
            if (object@paired) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for paired samples")
                int <- is_integer(object@ytrans[, 1L], ...)
                if (int)
                    SR_shift_1sample(object, fact = attr(int, "fact"))
                else
                    stop("cannot compute exact distribution with real-valued scores")
            } else if (is_2sample(object)) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for quadratic tests")
                else {
                    int <- is_integer(object@ytrans[, 1L], ...)
                    if (int)
                        SR_shift_2sample(object, fact = attr(int, "fact"))
                    else if (algorithm == "auto")
                        stop("split-up algorithm not implemented for quadratic tests")
                    else
                        stop("cannot compute exact distribution with real-valued scores")
                }
            } else
                stop(sQuote("object"), " is not a two-sample problem")
        }
)


### generic method for approximate null distributions
setGeneric("ApproxNullDistribution",
    function(object, ...) {
        standardGeneric("ApproxNullDistribution")
    }
)

### method for scalar test statistics
setMethod("ApproxNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        pls <- plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        ## <FIXME> can transform p, q, x instead of those </FIXME>
        pls <- sort((pls - expectation(object)) / sqrt(variance(object)))

        d <- function(x) {
            tmp <- abs(pls - x)
            mean(tmp == tmp[which.min(tmp)] & tmp < eps())
        }
        pvalue <- function(q) {
            switch(object@alternative,
                "less"      = mean(LE(pls, q)),
                "greater"   = mean(GE(pls, q)),
                "two.sided" = mean(GE(abs(pls), abs(q))))
        }
        pvalueinterval <- function(q, z = c(1, 0)) {
            pp <- if (object@alternative == "two.sided")
                      d(-q) + d(q) # both tails
                  else
                      d(q)
            pvalue(q) - z * pp
        }

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                mean(LE(pls, q))
            },
            q = function(p) {
                quantile(pls, probs = p, names = FALSE, type = 1L)
            },
            d = d,
            pvalue = function(q) {
                pvalue <- pvalue(q)
                attr(pvalue, "conf.int") <- confint_binom(round(pvalue * B), B)
                class(pvalue) <- "MCp"
                pvalue
            },
            midpvalue = function(q) {
                midpvalue <- pvalueinterval(q, z = 0.5)
                attr(midpvalue, "conf.int") <- confint_midp(round(midpvalue * B), B)
                class(midpvalue) <- "MCp"
                midpvalue
            },
            pvalueinterval = function(q) {
                setNames(pvalueinterval(q), nm = c("p_0", "p_1"))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else
                    sort(unique(drop(pls)))
            },
            name = "Monte Carlo Distribution")
    }
)

### method for max-type test statistics
setMethod("ApproxNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        pls <- plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        pls <- (pls - expect) / dcov

        ## <FIXME>
        ## pls is a rather large object (potentially)
        ## try not to copy it too often -- abs() kills you
        ## </FIXME>

        pmaxmin <- function() {
            pls <- switch(object@alternative,
                       "less"      = do.call("pmin.int", as.data.frame(t(pls))),
                       "greater"   = do.call("pmax.int", as.data.frame(t(pls))),
                       "two.sided" = do.call("pmax.int", as.data.frame(t(abs(pls)))))
            sort(pls)
        }

        d <- function(x) {
            tmp <- abs(pmaxmin() - x)
            mean(tmp == tmp[which.min(tmp)] & tmp < eps())
        }
        pvalue <- function(q) {
            switch(object@alternative,
                "less"      = mean(colSums(LE(pls, q)) > 0),
                "greater"   = mean(colSums(GE(pls, q)) > 0),
                "two.sided" = mean(colSums(GE(abs(pls), q)) > 0))
        }
        pvalueinterval <- function(q, z = c(1, 0)) {
            pp <- if (object@alternative == "two.sided")
                      d(-q) + d(q) # both tails
                  else
                      d(q)
            pvalue(q) - z * pp
        }

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                switch(object@alternative,
                    "less"      = mean(colSums(GE(pls, q)) == nrow(pls)),
                    "greater"   = mean(colSums(LE(pls, q)) == nrow(pls)),
                    "two.sided" = mean(colSums(LE(abs(pls), q)) == nrow(pls)))
            },
            q = function(p) {
                quantile(pmaxmin(), probs = p, names = FALSE, type = 1L)
            },
            d = d,
            pvalue = function(q) {
                pvalue <- pvalue(q)
                attr(pvalue, "conf.int") <- confint_binom(round(pvalue * B), B)
                class(pvalue) <- "MCp"
                pvalue
            },
            midpvalue = function(q) {
                midpvalue <- pvalueinterval(q, z = 0.5)
                attr(midpvalue, "conf.int") <- confint_midp(round(midpvalue * B), B)
                class(midpvalue) <- "MCp"
                midpvalue
            },
            pvalueinterval = function(q) {
                setNames(pvalueinterval(q), nm = c("p_0", "p_1"))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else
                    sort(unique(drop(pmaxmin())))
            },
            name = "Monte Carlo Distribution")
    }
)

### method for quad-type test statistics
setMethod("ApproxNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        pls <- plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        dcov <- object@covarianceplus
        expect <- expectation(object)
        a <- pls - expect
        pls <- sort(rowSums(crossprod(a, dcov) * t(a)))

        d <- function(x) {
            tmp <- abs(pls - x)
            mean(tmp == tmp[which.min(tmp)] & tmp < eps())
        }
        pvalue <- function(q) mean(GE(pls, q))
        pvalueinterval <- function(q, z = c(1, 0)) pvalue(q) - z * d(q)

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                mean(LE(pls, q))
            },
            q = function(p) {
                quantile(pls, probs = p, names = FALSE, type = 1L)
            },
            d = d,
            pvalue = function(q) {
                pvalue <- pvalue(q)
                attr(pvalue, "conf.int") <- confint_binom(round(pvalue * B), B)
                class(pvalue) <- "MCp"
                pvalue
            },
            midpvalue = function(q) {
                midpvalue <- pvalueinterval(q, z = 0.5)
                attr(midpvalue, "conf.int") <- confint_midp(round(midpvalue * B), B)
                class(midpvalue) <- "MCp"
                midpvalue
            },
            pvalueinterval = function(q) {
                setNames(pvalueinterval(q), nm = c("p_0", "p_1"))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else
                    sort(unique(drop(pls)))
            },
            name = "Monte Carlo Distribution")
    }
)


### S3 method for extraction of confidence intervals
confint.ScalarIndependenceTestConfint <-
    function(object, parm, level = 0.95, ...) {
        x <- if ("level" %in% names(match.call()))
                 object@confint(level)
             else
                 object@confint(object@conf.level)
        class(x) <- "ci"
        x
}
