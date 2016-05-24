# mosaic plot for a glm object
#  Allow it to use residuals of any type computed by residuals()
#  or to pass residuals calculated by another function, e.g., rstandard(), rstudent()
#
#  Allow to apply to any model with discrete factors
# last modified: 3/6/2009 1:51PM
#  - fixed buggy version using ideas from vcd:::plot.loglm
#  - now use $data component when it is a table

## TODO: move to utility.R
#is.discrete.model <- function(model)
# 	all(attr(terms(model), "dataClasses")[-1] %in% c("factor", "ordered"))


`mosaic.glm` <-	function(x, formula = NULL,
                         panel=mosaic, type=c("observed", "expected"),
                         residuals=NULL,
                         residuals_type = c("pearson", "deviance", "rstandard"),
                         gp = shading_hcl, gp_args = list(), ...)
{

	#require(vcd)
	if (!inherits(x,"glm")) stop("mosaic.glm requires a glm object")

	df.residual <- x$df.residual
	observed <- x$data

        if (is.null(formula) && inherits(observed, "table"))
            formula <- reformulate(names(dimnames(observed)))

        if (is.null(formula)) {
            if (is.environment(observed)) observed <- model.frame(x)
            else {
                if (!is.null(x$call$subset))
                    observed <- subset(observed, eval(x$call$subset, observed))
                if (!is.null(x$na.action))
                    observed <- observed[-x$na.action,]
            }
            ## get all factors excluding response
            factors <- sapply(observed, inherits, "factor")
            resp <- as.character(x$formula[[2]])
            factors <- observed[setdiff(colnames(observed[factors]), resp)]
            ## drop unused levels
            for(nm in names(factors)) {
                f <- factors[[nm]]
                if(is.factor(f) &&
                   length(unique(f[!is.na(f)])) < length(levels(f)))
                    factors[[nm]] <- factors[[nm]][, drop = TRUE]
            }
            ok <- TRUE
            ## check cross-classifying
            if (ok <- isTRUE(all(table(factors) == 1))) {
              warning("no formula provided, assuming ",
                      deparse(formula(terms(~ . , data = factors))),
                      "\n", call. = FALSE)
            }
            if (!ok)
              stop("cannot identify indexing factors from ", substitute(x),
                   "$data - please provide formula", call. = FALSE)
        }
        else {
            if (length(formula) == 3) formula <- formula[-2]
            ## get indexing factors allowing for missing data, subset etc
            factors <- do.call("model.frame", list(formula = formula,
                                                   data = observed,
                                                   subset = x$call$subset,
                                                   na.action = na.pass,
                                                   drop.unused.levels = TRUE))
            ## following loop needed due to bug in model.frame.default (fixed for R 2.12)
            for(nm in names(factors)) {
                f <- factors[[nm]]
                if(is.factor(f) && length(unique(f[!is.na(f)])) < length(levels(f)))
                    factors[[nm]] <- factors[[nm]][, drop = TRUE]
            }
            if (!is.null(x$na.action))
                factors <- factors[-x$na.action,]
        }

        if (x$family$family == "poisson") {
            observed <- as.table(tapply(x$y, factors, sum))
            expected <- as.table(tapply(fitted(x), factors, sum))
        }
        else{
            observed <- as.table(tapply(x$prior.weights, factors, sum))
            expected <- as.table(tapply(x$prior.weights * x$weights, factors, sum))
        }
        ## replace any missing values with zero
        observed[is.na(observed)] <- 0 #else strucplot would do this
        expected[is.na(expected)] <- 0

	type <- match.arg(tolower(type), c("observed", "expected"))
	if (any(observed < 0, na.rm = TRUE))
            stop("requires a non-negative response vector")

        ## reshape the residuals to conform to the structure of data

        ## if max one residual per cell, use residuals_type
        if (max(table(factors)) == 1) {
            residuals_type <- match.arg(tolower(residuals_type),
                                        c("pearson", "deviance", "rstandard"))
            if (missing(residuals))
                residuals <- if (residuals_type=="rstandard") rstandard(x)
                else residuals(x, type=residuals_type)
            residuals <- as.table(tapply(residuals, factors, sum))
            df <- x$df.residual
        }
        ## for marginal views, use aggregated working residuals
        else {
            residuals <- meanResiduals(x, factors)
            residuals_type <- "working" #what is this used for?
            df <- attr(residuals, "df")
            if (df == 0) {
                warning("There are zero degrees of freedom ",
                        "for the test of normality")
                df <- NA
            }
        }
        ## replace any missing values with zero
        residuals[is.na(residuals)] <- 0

	gp <- if (inherits(gp, "grapcon_generator"))
            do.call("gp", c(list(observed, residuals, expected, df),
                            as.list(gp_args)))
        else gp

	panel(observed, residuals=residuals, expected=expected, type=type,
              residuals_type=residuals_type, gp=gp, ...)
}

## convenience functions for sieve and assoc plots
sieve.glm <-
		function (x, ...)
{
	mosaic(x, panel = sieve, ...)
}

assoc.glm <-
		function (x, ...)
{
	mosaic(x, panel = assoc, ...)
}

