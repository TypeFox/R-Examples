# earth.varmod.R: build variance models for estimating prediction intervals
#
# TODO Extend the coverage table (print.inconf.tab) to show percentages in
# lower and upper intervals, so the user can check for asymmetry of the
# residuals.
#
# TODO Add QQ plot for prediction intervals, a "PIQ" plot
#
# TODO Consider making the code automatically detect non-monotonicity and
# issuing a warning.  Probably only possible for univariate models.
#
# TODO Could maybe prevent "Error in numericDeriv" by internally passing
# an explicit derivative in the call to nls.

VARMOD.METHODS <- c("const", "power", "power0",
                    "lm", "rlm", "earth", "gam",
                    "x.lm", "x.rlm", "x.earth", "x.gam")

TRACE.VARMOD         <- .3
TRACE.VARMOD.DETAILS <- .31 # will also cause plotting

# varmod returns a "varmod" object.
# y is the observed response (it is a n x 1 matrix)

varmod <- function(parent,
    method, exponent, conv, clamp, minspan,
    trace, x, y, model.var, ...)
{
    UseMethod("varmod")
}
varmod.earth <- function(parent,
    method, exponent, conv, clamp, minspan,
    trace, x, y, model.var, ...)
{
    check.classname(parent, substitute(parent), "earth")
    varmod.internal(parent,
                    method, exponent, conv, clamp, minspan,
                    trace, x, y, model.var, ...)
}
varmod.default <- function(parent,
    method, exponent, conv, clamp, minspan,
    trace, x, y, model.var, ...)
{
    warning0("varmod.default: varmods are not supported for \"",
             class(parent)[1], "\" objects\nContinuing anyway")

    varmod.internal(parent,
                    method, exponent, conv, clamp, minspan,
                    trace, x, y, model.var, ...)
}
varmod.internal <- function(parent,
    method, exponent=1, conv=1, clamp=.1,
    minspan=-5,
    trace=0, parent.x=NULL, parent.y=NULL,
    model.var, ...)
{
    # The following constant "lambda" was an argument to earth but I removed
    # it and hardcoded it here for simplicity in the earth interface.
    # We use lambda to transform the squared residuals as follows:
    #     transformed.resids = squared.resids ^ (lambda / 2)
    # So with lambda=1, we transform to absolute residuals, and if
    # lambda=2, then there is no transform.  We call the transformed
    # residuals the abs.resids in the code (which is actually the correct
    # nomenclature only when lambda is 1).  See also get.resids.name.
    lambda <- 1

    # likewise, rmethod is hardcoded here instead of being an arg to earth
    rmethod <- "hc12" # TODO doesn't match documentation "Variance models in earth"

    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    check.lambda.arg(lambda)
    check.exponent.arg(exponent, method)
    check.conv.arg(conv)
    check.clamp.arg(clamp)
    stopifnot(!is.null(parent.x))
    stopifnot(is.matrix(parent.x))
    stopifnot(!is.null(parent.y))
    stopifnot(is.matrix(parent.y))
    if(NCOL(parent.y) != 1)
        stop0("variance models are not supported for multiple response models")

    if(trace >= TRACE.VARMOD) {
        printf(
"\nvarmod method=\"%s\" rmethod=\"%s\" lambda=%g exponent=%g conv=%g clamp=%g minspan=%g:\n",
            method, rmethod, lambda, exponent, conv, clamp, minspan)
        if(trace == TRACE.VARMOD.DETAILS) {
            oldpar <- par(no.readonly=TRUE)
            on.exit(par(oldpar))
            par(mfrow=c(2, 3), mar=c(3, 3, 3, 1), mgp=c(1.5, .5, 0))
        }
    }
    n <- nrow(parent.x)
    df <- length(parent$selected.terms)
    leverages <- parent$leverages
    stopifnot(!is.null(leverages))
    leverages[leverages > .9] <- .9 # prevent any residual from being too influential
    correction <- switch(match.choices(rmethod,
                         c("hc0", "hc1", "hc2", "hc3", "hc12"), "varmod.rmethod"),
        hc0  = 1,
        hc1  = n / (n - df),
        hc2  = 1 / (1 - leverages),
        hc3  = 1 / (1 - leverages)^2,
        hc12 = n / ((n - df) * (1 - leverages)))

    squared.resids <- correction * (parent.y - predict(parent))^2 + model.var
    abs.resids <- squared.resids ^ (lambda / 2) # by default lambda=1, so this takes sqrt
    temp <- iterate.residmod(parent, abs.resids,
                             method, exponent, lambda, conv, clamp, minspan,
                             trace, parent.x, parent.y, ...)
        residmod  <- temp$residmod
        converged <- temp$converged
        iters     <- temp$iters

    if(trace >= TRACE.VARMOD)
        printf("\n")

    varmod             <- NULL
    varmod$call        <- make.call.generic(match.call(), "varmod")
    varmod$parent      <- parent
    varmod$method      <- method
    varmod$package     <- which.package(method)
    varmod$exponent    <- exponent
    varmod$lambda      <- lambda
    varmod$rmethod     <- rmethod
    varmod$converged   <- temp$converged
    varmod$iters       <- temp$iters
    varmod$residmod    <- residmod
    varmod$min.sd      <- get.min.sd(residmod, lambda, clamp)
    varmod$model.var   <- model.var
    varmod$abs.resids  <- abs.resids # transformed residuals (actually only abs when lambda is 1)
    varmod$parent.x    <- parent.x
    varmod$parent.y    <- parent.y
    class(varmod)      <- "varmod"
    varmod$iter.rsq    <- get.iter.rsq(varmod, abs.resids)
    varmod$iter.stderr <- get.iter.stderr(varmod, trace)
    attr(varmod, ".Environment") <- get.model.env(residmod, "varmod", trace)
    varmod
}
get.iter.rsq <- function(object, abs.resids) # return NULL if can't get rsq
{
    check.classname(object, substitute(object), "varmod")
    if(object$method == "const")
        return(NULL)
    fitted  <- fitted(object$residmod)
    weights <- weights(object$residmod)
    if(is.null(fitted) || is.null(weights))
        return(NULL)
    get.weighted.rsq(abs.resids, fitted, weights)
}
get.iter.stderr <- function(object, trace) # return if NULL if can't get stderr
{
    check.classname(object, substitute(object), "varmod")
    residmod <- object$residmod
    if(class(residmod)[1] %in% c("lm", "rlm", "nls")) {
        coef <- summary(residmod)$coefficients
        stopifnot(!is.null(coef[,"Std. Error"]))
        coef[,"Std. Error"]
    } else if(class(residmod)[1] == "gam" && object$package == "gam") {
        coef <- coefficients(summary.glm(residmod))
        stopifnot(!is.null(coef[,"Std. Error"]))
        coef[,"Std. Error"]
    } else if(class(residmod)[1] == "gam" && object$package == "mgcv") {
        # only the stderr for the intercept is available
        len.coef <- length(coefficients(residmod))
        std.err <- repl(NA, len.coef)
        std.err[1] <- summary(residmod)$p.table[1, "Std. Error"] # se of intercept
        std.err
    } else if(class(residmod)[1] == "earth") {
        trace2(trace, "--get.iter.stderr\n")
        y <- plotmo::plotmo_y(residmod, nresponse=1, trace)$y
        bx <- model.matrix(residmod)
        coef <- summary(lm(y~bx))$coefficients
        coef[,"Std. Error"]
    } else
        NULL
}
# iteratively reweighted least squares

issued.singularities.warning.global <- FALSE

iterate.residmod <- function(parent, abs.resids,
                            method, exponent, lambda, conv, clamp, minspan,
                            trace, parent.x, parent.y, ...)
{
    varmod   <- NULL
    max.iter <- 50
    weights  <- rep(1, nrow(parent.y))
    residmod <- NULL

    # following is needed because we want to issue singular warning at most once
    assignInMyNamespace("issued.singularities.warning.global", FALSE)

    # we always build the trace tab but only print it
    # if tracing is enabled or convergence failed
    trace.tab <- NULL

    for(iter in 1:max.iter) {
        residmod <- get.residmod(method, exponent, minspan, parent.x, parent.y,
                                 abs.resids, weights, trace, iter, parent, residmod, ...)

        coef.change <- get.coef.change(method, iter, residmod)

        # fill in enough of varmod for predict.varmod in get.residmod.weights
        varmod$parent   <- parent
        varmod$method   <- method
        varmod$exponent <- exponent
        varmod$lambda   <- lambda
        varmod$residmod <- residmod
        varmod$min.sd   <- get.min.sd(residmod, lambda, clamp)
        class(varmod)   <- "varmod"
        weights         <- get.residmod.weights(varmod, iter, trace)

        trace.tab <- update.trace.tab(trace.tab, trace, iter, max.iter,
                                      residmod, varmod, coef, coef.change,
                                      parent.y, weights, exponent)

        if(residmod.converged(coef.change, conv, iter, max.iter, method))
            break
    }
    converged <- residmod.converged(coef.change, conv, iter, max.iter, method)

    if(trace >= TRACE.VARMOD) {
         # || anyNA(coef(residmod)) || (conv >= 0 && !converged))
        if(trace == TRACE.VARMOD.DETAILS && inherits(residmod, "nls"))
            printf("\n")
        print(trace.tab[1:iter,], row.names=FALSE, digits=2)
    }
    if(!converged)
        warnf("varmod did not converge after %d iters (%s), final coefchange %.1f%%",
              iter, get.non.convergence.reason(trace.tab), mean(abs(coef.change)))

    list(residmod=residmod, converged=converged, iters=iter)
}
blank.plot <- function(main=NULL)
{
    plot(0, 0, col=0, bty="n", xlab="", ylab="",
         xaxt="n", yaxt="n", main=main)
}
trace.ncoef.global <- 0

update.trace.tab <- function(trace.tab, trace, iter, max.iter, residmod, varmod,
                             coef, coef.change, parent.y, weights, exponent)
{
    if(trace == TRACE.VARMOD.DETAILS) {
        if(inherits(residmod, "nls"))
            printf("\n")
        plotmo::plotmo(varmod, type="abs.residual",
            do.par=FALSE, degree1=1, degree2=0, trace=-1,
            pt.col=1, degree1.col=2, degree1.lwd=3,
            main="residmod first predictor")
    }
    coef <- coef(residmod)
    if(iter == 1) {
        # following needed because nbr of earth coefs can change
        assignInMyNamespace("trace.ncoef.global", length(coef))
        trace.tab <- as.data.frame(matrix(NA, nrow=max.iter, ncol=3+trace.ncoef.global))
        colnames(trace.tab) <- c("    iter", "weight.ratio", "coefchange%",
                 fix.coef.names(names(coef), colnames(parent.y), exponent))
    }
    trace.tab[iter,] <-
        c(iter, max(weights) / min(weights),
          mean(abs(coef.change)),
          c(coef, repl(NA, trace.ncoef.global))[1:trace.ncoef.global])
    trace.tab
}
# For debugging non-convergence.  In simulation about .4% of runs
# did not converge, mostly "oscillating", nearly all with sample
# sizes of less than 100.

get.non.convergence.reason <- function(trace.tab)
{

    nrow <- nrow(trace.tab)
    if(nrow < 7)
        return <- "short tab" # should never return this

    coefchange <- trace.tab[,"coefchange%"]
    c6 <- coefchange[nrow-6]
    c5 <- coefchange[nrow-5]
    c4 <- coefchange[nrow-4]
    c3 <- coefchange[nrow-3]
    c2 <- coefchange[nrow-2]
    c1 <- coefchange[nrow-1]
    c0 <- coefchange[nrow-0]

    # Example for oscillating:
    #    iter weight.ratio     coefchange% (Intercept)     y
    #      27           57   big   c3  3.7     0.00972 0.089
    #      28           63   small c2  3.4     0.00921 0.091
    #      29           57   big   c1  3.5     0.00971 0.089
    #      30           62   small c0  3.3     0.00922 0.091

    if(c0 < c1 && c1 > c2 && c2 < c3 && c3 > c0)
        return("oscillating-lo") # oscillating, last iter lower than prev

    else if(c0 > c1 && c1 < c2 && c2 > c3 && c3 < c0)
        return("oscillating-hi") # oscillating, last iter higher than prev

    if(c2 > c1 && c1 > c0) {        # only last two converged
        reason <- "converging2"
        if(c3 > c2) {               # only last three converged
            reason <- "converging3"
            if(c4 > c3) {           # only last four converged
                reason <- "converging4"
                if(c5 > c4) {       # only last five converged
                    reason <- "converging5"
                    if(c6 > c5)     # last six converged
                        reason <- "converging6"
                }
            }
        }
        return(reason)
    }
    if(c2 < c1 && c1 < c0) {        # last two diverged
        reason <- "diverging2"
        if(c3 < c2) {               # last three diverged
            reason <- "diverging3"
            if(c4 < c3) {           # last four diverged
                reason <- "diverging4"
                if(c5 < c4) {       # last five diverged
                    reason <- "diverging5"
                    if(c6 < c5)     # last six diverged
                        reason <- "diverging6"
                }
            }
        }
        return(reason)
    }
    "non-monotonic" # can find no other pattern
}
get.residmod <- function(method, exponent, minspan, parent.x, parent.y,
                         abs.resids, weights, trace, iter, parent, prev.residmod, ...)
{
    switch(method,
        const   = residmod.const(parent.x, parent.y, abs.resids,
                                 weights, trace),

        power   = residmod.power(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, parent, prev.residmod, iter, ...),

        power0  = residmod.power0(exponent, parent.x, parent.y, abs.resids,
                                  weights, trace, parent, prev.residmod, iter, ...),

        lm      = residmod.lm(exponent, parent.x, parent.y, abs.resids,
                              weights, trace, parent, ...),

        rlm     = residmod.rlm(exponent, parent.x, parent.y, abs.resids,
                               weights, trace, parent, ...),

        earth   = residmod.earth(exponent, minspan, parent.x, parent.y, abs.resids,
                                 weights, trace, parent, ...),

        gam     = residmod.gam(exponent, parent.x, parent.y, abs.resids,
                               weights, trace, parent, iter, ...),

        x.lm    = residmod.x.lm(exponent, parent.x, parent.y, abs.resids,
                                weights, trace, ...),

        x.rlm   = residmod.x.rlm(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, ...),

        x.earth = residmod.x.earth(exponent, minspan, parent.x, parent.y, abs.resids,
                                   weights, trace, ...),

        x.gam   = residmod.x.gam(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, iter, ...),

        stop0("illegal varmod.method \"", method, "\""))
}
prev.coef.global <- NULL

get.coef.change <- function(method, iter, residmod) # returns percents for each coef
{
    coef <- coef(residmod)

    # lm sometimes returns 2nd coef as NA if resids are all the same
    # TODO this should be an error? --- but that halts simulation tests
    if(anyNA(coef)) {
        if(!issued.singularities.warning.global) {
            warning0("singularities in residual model (coefs ",
                     paste.collapse(coef), ")")
            assignInMyNamespace("issued.singularities.warning.global", TRUE)
        }
        coef[is.na(coef)] <- 0
    }
    if(iter == 1)
        assignInMyNamespace("prev.coef.global", coef)

    if(method %in% c("earth", "x.earth")) # see comments in residmod.converged
        if(length(prev.coef.global) != length(coef))
            return(9999)

    coef.change <- abs(coef - prev.coef.global)
    prev <- abs(prev.coef.global)

    assignInMyNamespace("prev.coef.global", coef)

    # Divide by absolute value of previous coefficients.
    # But ensure no divide by near zero, by downweighting extremely small
    # coefs, thus preventing them from completely dominating the mean
    # change if the rest are large.
    # The 1e-8 prevents noise floor coefficients from preventing convergence.

    min <- max(.01 * max(prev), 1e-8)
    prev[prev < min] <- min
    100 * coef.change / prev # a percentage for each coef
}
residmod.converged <- function(coef.change, conv, iter, max.iter, method)
{
    if(conv < 0)
        iter >= abs(conv)
    else {
        method == "const" ||

        # TODO Since the earth basis funcs can change, looking at changes
        #      in the coefs can't be used to determine convergence.
        #      So for now, always do 1 iter for earth residual models.

        (method %in% c("earth", "x.earth")) ||

        # TODO following will sometimes create an intercept only model so unused
        # (method %in% c("earth", "x.earth") && iter >= 2) ||

        iter > 1 && mean(abs(coef.change)) < conv
    }
}
draw.residmod.weights <- function(w, main, min=NA, max=NA, median=NA) # for debugging
{
    plot(w, type="l", main=main, ylim=c(0, max(w, if(is.na(max)) 0 else max)))
    if(!is.na(min)) {
        abline(h=min, col=2)
        abline(h=max, col=2)
        abline(h=median, col=2)
    }
    else
        legend("topright", sprintf("max/min %.0f", max(w) / min(w)))
    lines(w) # replot over other annotations
}
# clamp to prevent extreme weights after squaring and inverse in get.residmod.weights

clamp.se <- function(se, iter, trace)
{
    median <- median(se)
    max.ratio <- 5 # 5 seems ok with (limited) simulation studies
    min <- median / max.ratio
    max <- max.ratio * median

    if(trace == TRACE.VARMOD.DETAILS)
        draw.residmod.weights(se, main=sprintf("iter %d: se", iter), min, max, median)

    se[se < min] <- min
    se[se > max] <- max

    se
}
# The variance for a regression on absolute residuals is proportional to
# the square of the regression model predicted value (Carrol and Ruppert
# book Section 3.3.3 and Table 3.3).

get.residmod.weights <- function(object, iter=0, trace=0)
{
    check.classname(object, substitute(object), "varmod")

    # square to convert se to variance, inverse to convert variance to weight
    weights <- 1 / clamp.se(predict.varmod(object, type="se"), iter, trace)^2

    # normalization is not strictly necessary, may help numerics
    weights <- weights / mean(weights)

    if(trace == TRACE.VARMOD.DETAILS)
        draw.residmod.weights(weights, "weights")

    weights
}
# We calculate lamba.factor.global only when necessary because the
# calculation can be slow.  Hence we need the following global variables.

lamba.global <- lamba.factor.global <- -999

update.lambda.factor <- function(lambda, trace)
{
    approx.equal <- function(x, y)
    {
        # allow for limited precision in doubles, also allows .33 for 1/3
        abs(x - y) < 1e-2
    }
    #--- update.lambda.factor starts here ---
    if(lambda != lamba.global) {
        assignInMyNamespace("lamba.global", lambda)

        # some values have been precalculated
        if(approx.equal(lambda, 2))
            assignInMyNamespace("lamba.factor.global", 1)
        # sqrt(pi / 2) = 1.2533, ratio mean dev to stddev, Geary 1935
        else if(approx.equal(lambda, 1))
            assignInMyNamespace("lamba.factor.global", sqrt(pi / 2))
        # (residuals^2)^(1/3) is approx normal by the Wilson-Hilferty
        # transform, although the left tail will still be short
        else if(approx.equal(lambda, 2/3))
            assignInMyNamespace("lamba.factor.global", 1.2464)
        else
        {
            rnorm(1) # seems to be necessary to make .Random.seed available
            old.seed <- .Random.seed
            set.seed(1) # for reproducibility
            # 1e6 below could be bigger but then slow
            assignInMyNamespace("lamba.factor.global",
                                1 / mean(rnorm(1e6)^2 ^ (lambda/2)))
            set.seed(old.seed)
        }
        if(trace >= TRACE.VARMOD)
            printf("lambda %g lamba.factor %g\n", lambda, lamba.factor.global)
    }
}
# scale a prediction by the residmod back to a standard deviation

to.sd <- function(abs.resids, lambda, trace=0)
{
    update.lambda.factor(lambda, trace)
    # pmax is necessary to prevent e.g. sqrt of neg prediction from residmod
    (lamba.factor.global * pmax(abs.resids, 0)) ^ (1 / lambda)
}
get.min.sd <- function(residmod, lambda, clamp=.1)
{
    predict <- predict(residmod)
    predict <- predict[predict > 0]
    stopifnot(length(predict) > 0)
    stopifnot(clamp >= 0, clamp <= 1)
    clamp * mean(to.sd(predict, lambda, 0))
}
check.lambda.arg <- function(lambda)
{
    check.numeric.scalar(lambda)
    if(lambda < .25 || lambda > 2)
        stop0("lambda=", lambda, " but it should be between 0.25 and 2")
}
# TRUE if estimation of variance depends only on the fitted response (not on x)
method.uses.fitted.response <- function(method)
{
    method %in% c("power", "power0", "lm", "rlm", "earth", "gam")
}
check.exponent.arg <- function(exponent, method)
{
    check.numeric.scalar(exponent)
    # TODO following restriction could be lifted but currently only partially implemented
    if(exponent != 1 && !method.uses.fitted.response(method))
        stop0("varmod.exponent argument is not allowed with method=\"", method, "\"\n",
              "(varmod.exponent is only allowed for varmod.methods that depend only ",
              "on the fitted response)")
    if(exponent < .1 || exponent > 5)
        stop0("varmod.exponent=", exponent, " but it should be between .1 and 5")
}
check.conv.arg <- function(conv)
{
    err <- function(conv)
        stop0("varmod.conv=", conv,
              " but it should be a negative integer ",
              "or a percent between 0 and 100")

    check.numeric.scalar(conv)
    if(conv < 0) {
        if(floor(conv) != conv) # conv is negative, check that it is an integer
            err(conv)
    } else if(conv == 0 || conv > 100)
        err(conv)
}
check.clamp.arg <- function(clamp)
{
    check.numeric.scalar(clamp)
    if(clamp < 0 || clamp > 1)
        stop0("varmod.clamp=", clamp, " but it should be between 0 and 1")
}
residmod.const <- function(parent.x, parent.y, abs.resids, weights, trace)
{
    # Predictions can be handled in a simple consistent way in
    # residmod.predict if instead of calculating the variance directly
    # here, we achieve the same result by building an intercept-only model
    # which always predicts mean(abs.resids).
    #
    # The conversion to a dataframe is necessary if the user later calls
    # plot(parent$varmod$residmod) or plotmo(parent$varmod$residmod).
    # Note that plotmo will call predict.varmod via predict.earth.

    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    lm(abs.resids~1, data=data, weights=weights, y=TRUE)
}
apply.exponent <- function(yhat, exponent, yname=colnames(yhat))
{
    stopifnot(!is.null(yname))
    check.vec(yhat, "yhat")
    # exponents of neg numbers are allowed only for integer exponents
    if(floor(exponent) != exponent) {
        check.that.most.are.positive(
            yhat, yname, sprintf("exponent=%g", exponent), "nonpositive")
        yhat[yhat < 0] <- 0 # don't want to take say sqrt of a neg number
    }
    yhat ^ exponent
}
nls.wrapper <- function(form, data, start, weights, abs.resids, trace)
{
    # We use algorithm="port" below because the default algorithm more often causes
    # "Error in numericDeriv: Missing value or an infinity produced"
    # Also, on test data we sometimes need more iterations than the default 50
    mod <- nls(formula=form,
               data=data, start=start, weights=weights,
               trace=(trace == TRACE.VARMOD.DETAILS),
               algorithm="port", control=list(maxiter=100))

    # make model data available for plotmo and plotres
    mod$x <- data[,-1,drop=FALSE]
    mod$y <- abs.resids

    # nls doesn't save the terms, so call$formula can confuse plotmo and plotres
    mod$call <- NULL

    mod
}
estimate.power.start.values <- function(prev.residmod, abs.resids, data, weights, trace, iter)
{
    if(is.null(prev.residmod)) { # first iteration in iterate.residmod?
        # use a linear model to estimate the start values
        lm <- lm(abs.resids~., data=data, weights=weights)
        coefs <- coef(lm)
        if(trace == TRACE.VARMOD.DETAILS) {
            plotmo::plotmo(lm, pt.col=2, do.par=FALSE, trace=-1,
                main=sprintf("iter 1: lm for start vals\nvarmod.method=\"power\""))
            plot(lm, which=1)
            blank.plot()
        }
        start <- list(coefs[1], coefs[2], exponent=1)
        if(trace >= TRACE.VARMOD)
            printf(
                "\n     start: (Intercept)=%.2g    coef=%.2g    exponent=%.2g\n\n",
                start[[1]], start[[2]], start[[3]])
    } else { # not first iteration
        # use previous model values as starting values
        coefs <- coef(prev.residmod)
        stopifnot(length(coefs) == 3)
        start <- list(coefs[1], coefs[2], coefs[3])
    }
    names(start) <- c("(Intercept)", "coef", "exponent")
    if(trace == TRACE.VARMOD.DETAILS)
        cat(sprintf("iter %d  RSS:   ", iter), names(start), "\n")
    start
}
residmod.power <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, parent, prev.residmod, iter, ...)
{
    if(exponent != 1) # TODO allow this?
        stop0("the exponent argument is not allowed with varmod.method=\"power\"")
    parent.fit <- predict(parent)
    check.that.most.are.positive(
        parent.fit, "predict(parent)", "varmod.method=\"power\"", "nonpositive")
    parent.fit[parent.fit < 0] <- 0 # force negative values to zero
    form <- abs.resids~`(Intercept)` + coef * RHS^exponent
    data <- data.frame(abs.resids, apply.exponent(parent.fit, exponent))
    colnames(data) <- c("abs.resids", "RHS")
    start <- estimate.power.start.values(prev.residmod, abs.resids,
                                         data, weights, trace, iter)
    nls.wrapper(form, data, start, weights, abs.resids, trace)
}
estimate.power0.start.values <- function(prev.residmod, abs.resids,
                                         data, weights, trace, iter)
{
    if(is.null(prev.residmod)) { # first iteration in iterate.residmod?
        # use a linear model to estimate the start values
        lm <- lm(abs.resids~.-1, data=data, weights=weights)
        coefs <- coef(lm)
        if(trace == TRACE.VARMOD.DETAILS) {
            plotmo::plotmo(lm, pt.col=2, do.par=FALSE, trace=-1,
                main=sprintf("iter 1: lm for start vals\nvarmod.method=\"power0\""))
            plot(lm, which=1)
            blank.plot()
        }
        start <- list(coefs[1], exponent=1)
        if(trace >= TRACE.VARMOD)
            printf(
                "\n     start: coef=%.2g    exponent=%.2g\n\n",
                start[[1]], start[[2]])
    } else { # not first iteration
        # use previous model values as starting values
        coefs <- coef(prev.residmod)
        stopifnot(length(coefs) == 2)
        start <- list(coefs[1], coefs[2])
    }
    names(start) <- c("coef", "exponent")
    if(trace == TRACE.VARMOD.DETAILS)
        cat(sprintf("iter %d  RSS:   ", iter), names(start), "\n")
    start
}
residmod.power0 <- function(exponent, parent.x, parent.y, abs.resids,
                            weights, trace, parent, prev.residmod, iter, ...)
{
    if(exponent != 1) # TODO allow this?
        stop0("the exponent argument is not allowed with varmod.method=\"power0\"")
    parent.fit <- predict(parent)
    check.that.most.are.positive(
        parent.fit, "predict(parent)", "varmod.method=\"power0\"", "nonpositive")
    parent.fit[parent.fit < 0] <- 0 # force negative values to zero
    data <- data.frame(abs.resids, apply.exponent(parent.fit, exponent))
    colnames(data) <- c("abs.resids", "RHS")
    start <- estimate.power0.start.values(prev.residmod, abs.resids,
                                          data, weights, trace, iter)
    nls.wrapper(abs.resids~coef * RHS^exponent,
                data, start, weights, abs.resids, trace)
}
residmod.lm <- function(exponent, parent.x, parent.y, abs.resids,
                        weights, trace, parent, ...)
{
    # we use RHS instead of colnames(parent.y) because we have applied exponent
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "RHS")
    lm(abs.resids~., data=data, weights=weights, y=TRUE)
}
residmod.rlm <- function(exponent, parent.x, parent.y, abs.resids,
                        weights, trace, parent, ...)
{
    # we use RHS instead of colnames(parent.y) because we have applied exponent
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "RHS")
    mod <- MASS::rlm(abs.resids~., data=data, weights=weights, method="MM")
    # make model data available for plotmo and plotres
    mod$data <- data
    mod
}
residmod.earth <- function(exponent, minspan, parent.x, parent.y, abs.resids,
                           weights, trace, parent, ...)
{
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "RHS")
    earth(abs.resids~., data=data, weights=weights,
          keepxy=TRUE, trace=trace, minspan=minspan, ...)
}
please.load.gam.package <- function()
{
    stop0("please load either the gam or mgcv package before using varmod.method=\"gam\"")
}
# Do we use the gam function in the gam or the mgcv package?
# Note that library(gam) has to be used before calling gam::gam, else the
# wrong "s" function is invoked.  This is because requireNamespace(gam)
# doesn't work there, even if we use gam::s when invoking gam.
# But CRAN check disallows library(gam) in the code (as from Jan 2015).
# So we have to ask the user to manually load the package if it is not loaded.

which.gam.package.is.loaded <- function()
{
    gam.package.loaded  <- "package:gam"  %in% search()
    mgcv.package.loaded <- "package:mgcv" %in% search()
    if(mgcv.package.loaded && gam.package.loaded) {
        # prevent downstream confusing error messages
        stop0("varmod.method=\"gam\" is not allowed when both the ",
               "'gam' and 'mgcv' packages are loaded")
    }
    if(gam.package.loaded)
        return("gam")
    if(mgcv.package.loaded)
        return("mgcv")
    please.load.gam.package()
}
which.package <- function(method)
{
    if(method %in% c("gam", "x.gam")) {
        if("package:gam" %in% search())
            return("gam")
        if("package:mgcv" %in% search())
            return("mgcv")
    }
    NULL
}
residmod.gam.aux <- function(form, data, weights, trace, iter)
{
    package.name <- which.gam.package.is.loaded()
    if(package.name == "gam") {
        if(trace >= TRACE.VARMOD && iter==1)
            printf("using the gam function from the 'gam' package\n")
        residmod <- gam::gam(formula=form, data=data, weights=weights)
        # We don't use x=TRUE else the x has colnames like s(x) which
        # confuses things later.  But we do save the data for plotmo.
        residmod$data <- data
    } else if(package.name == "mgcv") {
        if(trace >= TRACE.VARMOD && iter==1)
            printf("using the gam function from the 'mgcv' package\n")
        residmod <- mgcv::gam(formula=form, data=data, weights=weights)
        residmod$data <- data # for later access by plotmo etc.
    } else
        please.load.gam.package()
    residmod
}
residmod.gam <- function(exponent, parent.x, parent.y, abs.resids,
                         weights, trace, parent, iter, ...)
{
    form <- abs.resids ~ s(RHS)
    RHS <- apply.exponent(predict(parent), exponent)
    data <- data.frame(abs.resids, RHS)
    colnames(data) <- c("abs.resids", "RHS")
    residmod.gam.aux(form, data, weights, trace, iter)
}
residmod.x.lm <- function(exponent, parent.x, parent.y, abs.resids,
                          weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument is not allowed with varmod.method=\"x.lm\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    lm(abs.resids~., data=data, weights=weights, y=TRUE)
}
residmod.x.rlm <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument is not allowed with varmod.method=\"x.rlm\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    mod <- MASS::rlm(abs.resids~., data=data, weights=weights,
                     method="MM", y.ret=TRUE)
    # make model data available for plotmo and plotres
    mod$y <- abs.resids
    mod
}
residmod.x.earth <- function(exponent, minspan, parent.x, parent.y, abs.resids,
                             weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument is not allowed with varmod.method=\"x.earth\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    earth(abs.resids~., data=data, weights=weights,
          keepxy=TRUE, trace=trace, minspan=minspan, ...)
}
residmod.x.gam <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, iter, ...)
{
    if(exponent != 1)
        stop0("the exponent argument is not allowed with varmod.method=\"x.gam\"")
    if(ncol(parent.x) != 1)
        stop0("varmod.method=\"x.gam\" is not allowed when x has more than one column")
    form <- abs.resids ~ s(RHS)
    RHS <- parent.x[,1]
    data <- data.frame(abs.resids=abs.resids, RHS=RHS)
    colnames(data) <- c("abs.resids", "RHS")
    residmod.gam.aux(form, data, weights, trace, iter)
}
get.quant <- function(level) # e.g for level=.95 return 1.96
{
    check.level.arg(level, zero.ok=FALSE)
    stopifnot(level > 0, level < 1)
    level <- 1 - (1 - level) / 2 # .95 becomes .975
    qnorm(level)                 # .975 becomes 1.96
}
predict.se <- function(object, newdata)
{
    to.sd(predict.abs.residual(object, newdata), object$lambda)
}
predict.abs.residual <- function(object, newdata)
{
    # unfortunately needed to get model formulas to work for some varmod methods
    hack.colnames <- function(newdata, method)
    {
        if(!is.null(newdata) &&
           method %in% c("power", "power0", "lm", "rlm", "earth", "gam", "x.gam")) {
            if(NCOL(newdata) != 1) {
                stop0("predict.varmod: NCOL(newdata) must be 1 ",
                      "when method=\"", method, "\" (implementation restriction)")
            }
            newdata <- as.data.frame(newdata)
            colnames(newdata) <- "RHS"
        }
        newdata
    }
    if(is.null(newdata))
        abs.resid <- predict(object$residmod)
    else if(method.uses.fitted.response(object$method)) {
        parent.fit <- predict(object$parent, newdata=newdata)
        parent.fit <- apply.exponent(parent.fit, object$exponent)
        parent.fit <- data.frame(parent.fit)
        parent.fit <- hack.colnames(parent.fit, object$method)
        if(object$method %in% c("power", "power0"))
            parent.fit[parent.fit < 0] <- 0 # force negative values to zero
        abs.resid <- predict(object$residmod, newdata=parent.fit)
        stopifnot(length(abs.resid) == NROW(parent.fit))
    } else {
        newdata <- hack.colnames(newdata, object$method)
        abs.resid <- predict(object$residmod, newdata=newdata)
        stopifnot(length(abs.resid) == NROW(newdata))
    }
    abs.resid <- as.vector(abs.resid)

    # clamp at object$min.sd
    min.abs.resid <- (object$min.sd ^ object$lambda) / lamba.factor.global
    pmax(abs.resid, min.abs.resid)
}
get.parent.fit <- function(object, newdata)
{
    parent.fit <- predict(object$parent, newdata=newdata)
    check.vec(parent.fit, "parent.fit")
    stopifnot(!is.null(dim(parent.fit))) # check parent.fit is a matrix or dataframe
    parent.fit[,1]
}
predict.pint <- function(object, newdata, level) # newdata allowed
{
    se <- predict.se(object, newdata)
    parent.fit <- get.parent.fit(object, newdata)
    stopifnot(length(parent.fit) == length(se))
    quant <- get.quant(level)
    data.frame(fit = parent.fit,
               lwr = parent.fit - quant * se,
               upr = parent.fit + quant * se)
}
predict.cint <- function(object, newdata, level)
{
    if(!is.null(newdata))
        stop0("predict.varmod: newdata is not allowed with interval=\"cint\"")
    parent.fit <- get.parent.fit(object, newdata)
    se <- sqrt(object$model.var)
    stopifnot(length(se) == length(parent.fit))
    quant <- get.quant(level)
    data.frame(fit = parent.fit,
               lwr = parent.fit - quant * se,
               upr = parent.fit + quant * se)
}
predict.varmod <- function(
    object  = stop("no 'object' argument"),
    newdata = NULL,
    type    = c("pint", "cint", "se", "abs.residual"),
    level   = .95,
    trace   = FALSE, # unused but needed for plotmo
    ...)
{
    check.level95 <- function(level)
    {
        check.level.arg(level, zero.ok=TRUE)
        if(level != .95)
            stop0("predict.varmod: the level argument is not allowed with type=\"",
                  type, "\"")
    }
    check.classname(object, substitute(object), "varmod")
    warn.if.dots(...)
    switch(match.arg1(type, "type"),
        pint = predict.pint(object, newdata, level),
        cint = predict.cint(object, newdata, level),
        se   = {
            check.level95(level)
            predict.se(object, newdata)
        },
        abs.residual = {
            check.level95(level)
            predict.abs.residual(object, newdata)
        })
}
# Example: if digits=3, then "%.*f" becomes "%.3f"
# Needed because R printf doesn't support * in printf formats
# and we need it to make the digits arg work in printfs

dot.star.to.digits <- function(s, digits)
{
    check.integer.scalar(digits, min=1)
    stopifnot(floor(digits) == digits)
    stopifnot(digits > 0, digits < 20)
    gsub(".*", sprintf(".%d", digits), s, fixed=TRUE)
}
# Example:
#       fix.coef.names(coef.names=h(y-123), resp.name="y", exponent=.5)
#   returns
#       h(sqrt(y)-123) # the func knows that the special case of exponent=.5 is sqrt

fix.coef.names <- function(coef.names, resp.name, exponent)
{
    if(length(coef.names) == 1)
        return(coef.names) # do nothing if intercept only model
    stopifnot(length(resp.name) == 1)
    stopifnot(exponent > 0)
    new.resp.name <-
        if(exponent > .33 && exponent < .34)
            sprintf("cbrt(%s)", resp.name)
        else if(exponent == .5)
            sprintf("sqrt(%s)", resp.name)
        else if(exponent == 1)
            resp.name
        else if(exponent == 2)
            sprintf("sq(%s)", resp.name)
        else
            sprintf("%s^%.3g", resp.name, exponent)
    coef.names <- gsub("RHS", resp.name, coef.names, fixed=TRUE)
    if(exponent == 1)
        coef.names
    else {
        # TODO revisit, will fail if resp.name is substring of a token in
        #      coef.names or if resp.name="h" and coef.names="h(h-12)"
        gsub(resp.name, new.resp.name, coef.names, fixed=TRUE)
    }
}
# restore original exponent, it doesn't get scaled like the other coefficients

restore.exponent <- function(coef, org.coef, method)
{
    if(method == "power")
        coef[3] <- org.coef[3]  # exponent is in coef[3]
    else if(method == "power0")
        coef[2] <- org.coef[2]  # exponent is in coef[2]
    coef
}
coef.varmod <- function(object, as.sd=TRUE, ...)
{
    warn.if.dots(...)
    coef <- coef(object$residmod)
    if(is.null(coef))
        stop0("coef.varmod: cannot get coefficients for \"",
               class(object$residmod)[1], "\" residmod")
    as.sd <- check.boolean(as.sd)
    if(as.sd) {
        org.coef <- coef
        negs <- coef < 0
        coef <- to.sd(abs(coef), object$lambda)
        coef[negs] <- -coef[negs]
        coef <- restore.exponent(coef, org.coef, object$method)
    }
    names(coef) <- fix.coef.names(names(coef),
                                  colnames(object$parent.y), object$exponent)
    coef
}
VARMOD.COEF.TAB.STYLES <- c("standard", "unit")

get.varmod.coef.tab <- function(
    object,
    style  = VARMOD.COEF.TAB.STYLES)
{
    style <- match.arg1(style, "style")
    coef <- coef.varmod(object, as.sd=TRUE)

    # if style="unit", normalize coef if possible
    unit <- 1
    if(style == "unit") {
        # choose which coef will be the unit
        if(length(coef) == 1) # method == "const"?
            unit <- abs(coef[1])
        else
            unit <- abs(coef[2])
        if(unit < 1e-3) {
            warning0("coef=", unit, " is very small, forcing style=\"standard\"\n")
            style <- "standard"
            unit <- 1
        }
    }
    org.coef <- coef
    coef <- coef / unit
    coef <- restore.exponent(coef, org.coef, object$method)

    # get stderr
    NAs <- repl(NA, length(coef))
    stderr <- NAs
    if(!is.null(object$iter.stderr)) {
        org.stderr <- object$iter.stderr
        stderr <- to.sd(org.stderr, object$lambda) / unit
        stderr <- restore.exponent(stderr, org.stderr, object$method)
    }
    abs.coef <- abs(coef)
    stderr.percent <- 100 * stderr / abs.coef
    coef.names <- names(coef)
    coef.names <- gsub("`", "", coef.names) # remove backquotes added by lm etc.
    if(style == "unit") {
        coef.tab <- data.frame(c(coef, unit), c(stderr, NA), c(stderr.percent, NA))
        rownames(coef.tab) <- c(coef.names, "unit")
    } else {
        coef.tab <- data.frame(coef, stderr, stderr.percent)
        rownames(coef.tab) <- coef.names
    }
    if(object$method %in% c("earth", "x.earth")) {
        order <- reorder.earth(object$residmod, decomp="anova")
        coef.tab <- coef.tab[order,]
    }
    colnames(coef.tab) <- c("coefficients", "iter.stderr", "iter.stderr.percent")
    coef.tab
}
get.interval.tab <- function(object, level)
{
    level <- check.level.arg(level, zero.ok=FALSE)
    predict <- predict.varmod(object, type="pint", level=level)
    interval <- predict$upr - predict$lwr
    interval <- interval[order(interval)]
    tab <- data.frame(
        " ", mean(interval),
        " ", interval[1],
        " ", interval[length(interval)],
        " ", interval[length(interval)] / interval[1])
    colnames(tab) <- c(
        " ", "mean",
        " ", "smallest",
        " ", "largest",
        " ", "ratio")
    rownames(tab) <- sprintf("%g%% prediction interval", 100*level)
    tab
}
percent.inconf <- function(object, level, parent.y, newdata)
{
    predict <- predict.varmod(object, newdata, type="pint", level=level)
    inconf <- parent.y >= predict$lwr & parent.y <= predict$upr
    100 * sum(inconf) / length(inconf)
}
print.inconf.tab <- function(object, parent.y, newdata)
{
    if(NCOL(parent.y) != 1) {
        warning0("multiple response model: the table is for the first response")
        parent.y <- parent.y[,1]
    }
    stopifnot(is.numeric(parent.y) || is.logical(parent.y))

    inconf68 <- percent.inconf(object, .68, parent.y, newdata)
    inconf80 <- percent.inconf(object, .80, parent.y, newdata)
    inconf90 <- percent.inconf(object, .90, parent.y, newdata)
    inconf95 <- percent.inconf(object, .95, parent.y, newdata)

    # .5 below adjusts for rounding in printf %.0f
    lt <- function(x, level) if(x < level-.5) "<" else " "

    tab <- data.frame(
        " ", sprintf("%.0f%s ", inconf68, lt(inconf68, 68)),
        " ", sprintf("%.0f%s ", inconf80, lt(inconf80, 80)),
        " ", sprintf("%.0f%s ", inconf90, lt(inconf90, 90)),
        " ", sprintf("%.0f%s ", inconf95, lt(inconf95, 95)))

    colnames(tab) <- c(
        " ", "68% ",
        " ", "80% ",
        " ", "90% ",
        " ", "95% ")

    if(is.null(newdata))
        rowname <- "response values in prediction interval"
    else
        rowname <- "newdata in prediction interval"
    rownames(tab) <- rowname
    print(tab)

    # return value is the table but not in string form
    tab <- data.frame(inconf68, inconf80, inconf90, inconf95)
    colnames(tab) <- c("68%", "80%", "90%", "95%")
    rownames(tab) <- rowname
    tab
}
print.varmod <- function(
    x       = stop("no 'x' argument"), # x is a varmod object
    level   = .95,        # use 0 to not print the interval tabs
    style   = "standard", # one of VARMOD.COEF.TAB.STYLES
    digits  = 2,
    newdata = NULL,
    ...)
{
    check.classname(x, substitute(x), "varmod")
    object <- x # minimize confusion with x, the regression input matrix
    remove(x)   # not necessary but prevents mistakes later
    warn.if.dots(...)
    if(!is.null(newdata)) { # if newdata, print just the inconf table
        object$inconf.tab <-
            print.inconf.tab(object,
                plotmo::plotmo_response(object$parent, newdata, trace=0, ...),
                newdata)
        return(invisible(object))
    }
    printf("method \"%s\"", object$method)
    if(!is.null(object$package))
        printf(" (%s package)", object$package)
    space <- if(object$exponent != 1 || object$lambda != 1) "" else "  "
    if(object$exponent != 1)
        printf("%s  exponent %.3f", space, object$exponent)
    if(object$lambda != 1)
        printf("%s  lambda %g", space, object$lambda)
    printf("%s  min.sd %.3g", space, object$min.sd)
    if(!is.null(object$iter.rsq)) {
        printf("%s  iter.rsq %.3f", space, object$iter.rsq)
        # TODO prints too many digits
        # printf(dot.star.to.digits(", iter.rsq %.*f", digits+1), object$iter.rsq)
    }
    # coef tab
    printf("\n\nstddev of predictions%s:\n",
           if(style == "unit") " (scaled by unit)" else "")
    tab <- object$coef.tab
    if(is.null(tab)) { # needed if did not come here via summary.varmod
        tab <- get.varmod.coef.tab(object, style)
        object$coef.tab <- tab # for return value of this function
    }
    tab$coefficients <- zapsmall(tab$coefficients, digits+1)
    # sprintf below so print "NA" not "<NA>"
    tab$iter.stderr  <- sprintf("%g", zapsmall(tab$iter.stderr, digits))

    # convert iter.stderr.percent to character and print "big" if appropriate
    tab$iter.stderr.percent.as.char <- sprintf("%.0f", tab$iter.stderr.percent)
    tab$iter.stderr.percent.as.char[tab$iter.stderr.percent >= 1e3] <- "big"
    tab$iter.stderr.percent <- NULL
    colnames(tab) <- c("coefficients", "iter.stderr", "iter.stderr%")
    print(tab, digits=digits)

    # interval and inconf tabs
    level <- check.level.arg(level, zero.ok=TRUE)
    if(is.specified(level)) {
        stopifnot(level == object$level)
        printf("\n")
        tab <- object$interval.tab
        if(is.null(tab)) {
            tab <- get.interval.tab(object, level)
            object$interval.tab <- tab # for return value of this function
        }
        print(tab, digits=digits)
        printf("\n")
        object$inconf.tab <- print.inconf.tab(object, object$parent.y, newdata=NULL)
    }
    invisible(object)
}
print.summary.varmod <- function(
    x       = stop("no 'x' argument"), # x is a summary.varmod object
    level   = x$level,
    style   = x$style,
    digits  = x$digits,
    newdata = x$newdata,
    ...)
{
    check.classname(x, substitute(x), "varmod")
    warn.if.dots(...)
    if(is.null(level))
        level <- .95
    if(is.null(style))
        style <- "standard"
    if(is.null(digits))
        digits <- 2
    if(!is.null(newdata)) # if newdata, print just the inconf table
        print.varmod(x, level, style, digits, newdata)
    else {
        printcall("Parent model: ", x$parent$call)
        printf("\n")
        print.varmod(x, level, style, digits)
        printf("\nRegression submodel (%s):\n", get.resids.name(x))
        if(!(class(x$residmod)[1] %in% c("lm", "x.lm")))
            printf("\n")
        print(x$residmod, digits=digits)
    }
    invisible(x)
}
summary.varmod <- function(
    object  = stop("no 'object' argument"),
    level   = .95,
    style   = "standard", # one of VARMOD.COEF.TAB.STYLES
    digits  = 2,
    newdata = NULL,
    ...)
{
    check.classname(object, substitute(object), "varmod")
    warn.if.dots(...)
    object$level        <- level   # pass level on to print.summary.varmod
    object$style        <- style   # ditto
    object$digits       <- digits  # ditto
    object$newdata      <- newdata # ditto
    object$coef.tab     <- get.varmod.coef.tab(object, style)
    object$interval.tab <- get.interval.tab(object, level)
    # TODO add inconf table here too
    class(object)   <- c("summary.varmod", "varmod")
    object
}
get.resids.name <- function(object)
{
    if(object$lambda == 1)
        sprintf("Abs Residuals")
    else if(object$lambda == 2)
        sprintf("Squared Residuals")
    else
        sprintf("(Squared Residuals) ^ %.2g", object$lambda / 2)
}
get.varmod.ylab <- function(object, as.sd)
{
    sprintf("ParentMod %s", if(as.sd) "StdDev" else get.resids.name(object))
}
min.sd.line <- function(object, min.sd.col, lwd) # draw horizontal line at min.sd
{
    if(is.specified(min.sd.col)) {
        # TODO need to apply lambda exponent here?
        abline(h=object$min.sd / lamba.factor.global, col=min.sd.col, lty=2, lwd=lwd)
    }
}
sd.axis <- function(object) # draw righthand axis in standard deviation scale
{
    # for righthand axis
    sd <- lamba.factor.global * object$abs.resids ^ (1 / object$lambda)
    pretty.sd <- pretty(range(sd))
    axis(side=4, at=pretty.sd / lamba.factor.global, labels=pretty.sd, srt=90)
    # the line setting depends on the axis margin lines
    mtext(get.varmod.ylab(object, as.sd=TRUE), side=4, cex=par("cex"),
          line=if(par("mgp")[1] < 1.8) 1.4 else 1.8)
}
plot.varmod <- function(
    x          = stop("no 'x' argument"),
    which      = 1:4,
    do.par     = NULL,
    info       = FALSE,
    cex        = NULL,
    caption    = NULL,
    line.col   = 2,
    min.sd.col = line.col,
    trace      = 0,
    ...)    # unused, for compat with the generic
{
    check.classname(x, substitute(x), "varmod")
    object <- x # minimize confusion with x, the regression input matrix
    remove(x)   # needed else plotmo.x gets this x instead of the x matrix
    warn.if.dots(...)
    trace <- as.numeric(check.integer.scalar(trace, logical.ok=TRUE))
    info <- check.boolean(info)
    check.index(which, "which", 1:4)
    do.par <- check.do.par(do.par, length(which)) # do.par is 0, 1, or 2
    # prepare caption --- we need it now for do.par() but
    # can only display it later after at least one plot
    stopifnot.string(caption, allow.empty=TRUE, null.ok=TRUE)
    if(length(which) > 1 && do.par && is.null(caption)) # auto caption?
        caption <- sprintf("Variance Model  method=\"%s\"\nParentMod: %s",
                           object$method, strip.deparse(object$parent$call))
    main <- dot("main", ...)
    if(do.par) {
        oldpar <- par(no.readonly=TRUE)
        do.par(nfigs=length(which), caption=caption, main1=main,
               xlab1=NULL, ylab1=NULL, trace=trace,
               nlines.in.main=2, def.right.mar=3,
               def.font.main=1, # for compat with lm.plot
               ...)
        if(do.par == 1)
            on.exit(par(oldpar), add=TRUE)
    } else { # do.par=FALSE
        oldpar <- do.par.dots(..., trace=trace)
        if(length(oldpar))
            on.exit(do.call(par, oldpar), add=TRUE)
    }
    if(is.null(cex))
        cex <- pt.cex(length(object$parent.y))
    if(is.specified(main))
        main <- repl(main, 4) # recycle for up to 4 plots
    ylim <- fix.lim(c(min(object$abs.resids, 0), max(object$abs.resids)))
    parent.fitted <- predict(object$parent)[,1]
    order <- order(parent.fitted)
    smooth.col <- if(info) 2 else 0
    lwd <- 1
    for(iwhich in seq_along(which)) {
        if(which[iwhich] == 1) {            #--- fitted vs parent fitted ---
            plot(parent.fitted[order], object$abs.resids[order],
                 main=if(!is.specified(main))
                        sprintf("%s vs Fitted", get.resids.name(object))
                      else
                        main[iwhich],
                 ylim=ylim, pch=20, cex=cex, xlab="Fitted",
                 ylab=get.varmod.ylab(object, as.sd=FALSE))
            min.sd.line(object, min.sd.col, lwd) # horizontal line at min.sd
            sd.axis(object)                # right hand axis in stddev scale
            # fitted values of residual model
            fitted <- predict.varmod(object, type="abs.residual")
            lines(parent.fitted[order], fitted[order], col=line.col, lwd=lwd)
            if(info) {
                # lowess smooth
                smooth <- lowess(parent.fitted[order],
                                 object$abs.resids[order], f=.5)
                lines(smooth$x, smooth$y, col=smooth.col, lwd=1)
            }

        } else if(which[iwhich] == 2) {     #--- fitted vs parent first pred ---
            plotmo::plotmo(object, type="abs.residual",
                ylim=ylim, degree1=1, degree2=0, do.par=FALSE,
                trace=if(trace==0) -1 else trace,
                pt.col=1, pt.cex=cex, degree1.col=line.col,
                degree1.lwd=lwd, smooth.col=smooth.col,
                ylab=get.varmod.ylab(object, as.sd=FALSE),
                main=if(!is.specified(main))
                        sprintf("%s vs First Predictor", get.resids.name(object))
                     else
                        main[iwhich])
            min.sd.line(object, min.sd.col, lwd) # horizontal line at min.sd
            sd.axis(object)                 # right hand axis in stddev scale

        } else if(which[iwhich] == 3) {     #--- residual plot ---
            plotmo::plotres(object$residmod, which=3,
                            do.par=FALSE, center=FALSE,
                            xlab=get.varmod.ylab(object, as.sd=FALSE),
                            ylab="VarMod Residuals", info=info)

        } else if(which[iwhich] == 4) {     #--- model selection graph ---
            if(class(object$residmod)[1] == "earth")
                plot.earth(object$residmod, which=1, do.par=FALSE,
                           main=if(!is.specified(main)) "VarMod Model Selection"
                                else                    main[iwhich])
        } else
            stop0("plot.varmod: illegal value %g in 'which' argument",
                   which[iwhich])
    }
    draw.caption(caption, ...)
    invisible()
}
