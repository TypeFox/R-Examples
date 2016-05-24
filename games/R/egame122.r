##' @export
predict.egame122 <- function(object, newdata, type = c("outcome", "action"),
                             na.action = na.pass, ...)
{
    type <- match.arg(type)

    if (missing(newdata) || is.null(newdata)) {
        ## use original data if 'newdata' not supplied
        mf <- object$model
    } else {
        ## get rid of left-hand variables in the formula, since they're not
        ## needed for fitting
        formulas <- Formula(delete.response(terms(formula(object$formulas))))

        mf <- model.frame(formulas, data = newdata, na.action = na.action,
                          xlev = object$xlevels)

        ## check that variables are of the right classes
        Terms <- attr(object$model, "terms")
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, mf)
    }

    regr <- list()
    for (i in seq_len(length(object$formulas)[2]))
        regr[[i]] <- model.matrix(object$formulas, data = mf, rhs = i)

    ## get action probabilities, as given by fitted model parameters
    ans <- makeProbs122(object$coefficients, regr = regr, link = object$link, type
                        = object$type)

    if (type == "outcome") {
        ans <- data.frame(actionsToOutcomes122(ans, log.p = FALSE))
        names(ans) <- paste("Pr(", levels(object$y), ")", sep = "")
    } else {
        ans <- as.data.frame(ans)
        nl <- paste(levels(object$y)[1:2], collapse = " or ")
        nr <- paste(levels(object$y)[3:4], collapse = " or ")
        names(ans)[1:2] <- paste("Pr(", c(nl, nr), ")", sep = "")
        names(ans)[3:4] <- paste("Pr(", levels(object$y)[1:2], "|", nl, ")",
                                 sep="")
        names(ans)[5:6] <- paste("Pr(", levels(object$y)[3:4], "|", nr, ")",
                                 sep="")
    }

    return(ans)
}

sbi122 <- function(y, regr, link)
{
    names(regr) <- character(length(regr))
    names(regr)[1:6] <- c("X1", "X2", "X3", "X4", "Z2", "Z4")

    if (link == "probit") {
        fam <- binomial(link = "probit")
    } else {
        fam <- binomial(link = "logit")
    }

    ## regression for player 2's choice after 1 moves "left"
    ZL <- regr$Z2[y == 1 | y == 2, , drop = FALSE]
    yL <- as.numeric(y == 2)[y == 1 | y == 2]
    mL <- suppressWarnings(glm.fit(ZL, yL, family = fam))
    p2 <- as.numeric(regr$Z2 %*% coef(mL))
    p2 <- if (link == "probit") pnorm(p2) else plogis(p2)

    ## regression for player 2's choice after 1 moves "right"
    ZR <- regr$Z4[y == 3 | y == 4, , drop = FALSE]
    yR <- as.numeric(y == 4)[y == 3 | y == 4]
    mR <- suppressWarnings(glm.fit(ZR, yR, family = fam))
    p4 <- as.numeric(regr$Z4 %*% coef(mR))
    p4 <- if (link == "probit") pnorm(p4) else plogis(p4)

    ## regression for player 1's choice
    X1 <- cbind(-(1-p2) * regr$X1, -p2 * regr$X2, (1 - p4) * regr$X3, p4 *
                regr$X4)
    y1 <- as.numeric(y == 3 | y == 4)
    m1 <- glm.fit(X1, y1, family = fam)

    ## need to multiply by sqrt(2), see comments on 'sbi12' in 'egame12.r'
    ans <- sqrt(2) * c(coef(m1), coef(mL), coef(mR))
    return(ans)
}

makeSDs122 <- function(b, regr, type)
{
    sds <- vector("list", 8)
    rcols <- sapply(regr, ncol)

    if (length(rcols) == 7) {  ## sdByPlayer == FALSE
        v <- exp(as.numeric(regr[[7]] %*% b))
        for (i in 1:8) sds[[i]] <- v
    } else {
        v1 <- exp(as.numeric(regr[[7]] %*% b[1:rcols[7]]))
        v2 <- exp(as.numeric(regr[[8]] %*% b[(rcols[7]+1):length(b)]))
        if (type == "agent") {
            sds[[1]] <- sds[[2]] <- v1
            sds[[3]] <- sds[[4]] <- sds[[5]] <- sds[[6]] <- v2
        } else {
            sds[[1]] <- sds[[2]] <- sds[[3]] <- sds[[4]] <- v1
            sds[[5]] <- sds[[6]] <- sds[[7]] <- sds[[8]] <- v2
        }
    }

    return(sds)
}

makeProbs122 <- function(b, regr, link, type)
{
    private <- type == "private"

    utils <- makeUtils(b, regr, nutils = 6,
                       unames = c("u11", "u12", "u13", "u14", "u22", "u24"))

    ## length(utils$b) == 0 means no terms left for the variance components, so
    ## set these to 1
    if (length(utils$b) == 0) {
        sds <- as.list(rep(1, 6))
    } else {
        sds <- makeSDs122(utils$b, regr, type)
    }

    linkfcn <- switch(link,
                      logit = function(x, sd = 1) plogis(x, scale = sd),
                      probit = pnorm)

    if (private) {
        sd6 <- sqrt(sds[[7]]^2 + sds[[8]]^2)
    } else {
        sd6 <- sqrt(sds[[5]]^2 + sds[[6]]^2)
    }
    p6 <- finiteProbs(linkfcn(utils$u24, sd = sd6))
    p5 <- 1 - p6

    if (private) {
        sd4 <- sqrt(sds[[5]]^2 + sds[[6]]^2)
    } else {
        sd4 <- sqrt(sds[[3]]^2 + sds[[4]]^2)
    }
    p4 <- finiteProbs(linkfcn(utils$u22, sd = sd4))
    p3 <- 1 - p4

    sd2 <- if (private) {
        sqrt(p3^2 * sds[[1]]^2 + p4^2 * sds[[2]]^2 + p5^2 * sds[[3]]^2 +
             p6^2 * sds[[4]]^2)
    } else sqrt(sds[[1]]^2 + sds[[2]]^2)
    p2 <- p5 * utils$u13 + p6 * utils$u14 - p3 * utils$u11 - p4 * utils$u12
    p2 <- finiteProbs(linkfcn(p2, sd = sd2))
    p1 <- 1 - p2

    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}

actionsToOutcomes122 <- function(probs, log.p = TRUE)
{
    ans <- cbind(log(probs$p1) + log(probs$p3),
                 log(probs$p1) + log(probs$p4),
                 log(probs$p2) + log(probs$p5),
                 log(probs$p2) + log(probs$p6))
    if (!log.p) ans <- exp(ans)
    return(ans)
}

logLik122 <- function(b, y, regr, link, type, ...)
{
    names(regr) <- character(length(regr))
    names(regr)[1:6] <- c("X1", "X2", "X3", "X4", "Z2", "Z4")

    probs <- makeProbs122(b, regr, link, type)
    logProbs <- actionsToOutcomes122(probs, log.p = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad122 <- function(b, y, regr, link, type, ...)
{
    names(regr) <- character(length(regr))
    names(regr)[1:6] <- c("X1", "X2", "X3", "X4", "Z2", "Z4")

    u <- makeUtils(b, regr, nutils = 6,
                   unames = c("u11", "u12", "u13", "u14", "u22", "u24"))
    p <- makeProbs122(b, regr, link, type)
    rcols <- sapply(regr, ncol)
    n <- nrow(regr$X1)

    if (link == "probit" && type == "private") {
        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg2 <- dnorm(u$u22 / sqrt(2)) * (regr$Z2 / sqrt(2))
        dp4dg4 <- matrix(0L, nrow = n, ncol = rcols[6])
        dp4 <- cbind(dp4db, dp4dg2, dp4dg4)
        dp3 <- -dp4

        dp6db <- dp4db
        dp6dg2 <- matrix(0L, nrow = n, ncol = rcols[5])
        dp6dg4 <- dnorm(u$u24 / sqrt(2)) * (regr$Z4 / sqrt(2))
        dp6 <- cbind(dp6db, dp6dg2, dp6dg4)
        dp5 <- -dp6

        num2 <- p$p5 * u$u13 + p$p6 * u$u14 - p$p3 * u$u11 - p$p4 * u$u12
        denom2 <- sqrt(p$p3^2 + p$p4^2 + p$p5^2 + p$p6^2)
        dn2 <- dnorm(num2 / denom2)
        dp2db1 <- dn2 * (-p$p3) * regr$X1 / denom2
        dp2db2 <- dn2 * (-p$p4) * regr$X2 / denom2
        dp2db3 <- dn2 * p$p5 * regr$X3 / denom2
        dp2db4 <- dn2 * p$p6 * regr$X4 / denom2
        dp2dg2 <- dn2 * ((u$u11 - u$u12) * denom2 - num2 * (p$p4 - p$p3) / denom2)
        dp2dg2 <- (dp2dg2 * dp4dg2) / denom2^2
        dp2dg4 <- dn2 * ((u$u14 - u$u13) * denom2 - num2 * (p$p6 - p$p5) / denom2)
        dp2dg4 <- (dp2dg4 * dp6dg4) / denom2^2
        dp2 <- cbind(dp2db1, dp2db2, dp2db3, dp2db4, dp2dg2, dp2dg4)
        dp1 <- -dp2
    } else if (type == "agent") {
        dlink <- switch(link,
                        logit = dlogis,
                        probit = dnorm)

        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg2 <- dlink(u$u22 / sqrt(2)) * regr$Z2 / sqrt(2)
        dp4dg4 <- matrix(0L, nrow = n, ncol = rcols[6])
        dp4 <- cbind(dp4db, dp4dg2, dp4dg4)
        dp3 <- -dp4

        dp6db <- dp4db
        dp6dg2 <- matrix(0L, nrow = n, ncol = rcols[5])
        dp6dg4 <- dlink(u$u24 / sqrt(2)) * regr$Z4 / sqrt(2)
        dp6 <- cbind(dp6db, dp6dg2, dp6dg4)
        dp5 <- -dp6

        dn2 <- dlink((p$p5 * u$u13 + p$p6 * u$u14 - p$p3 * u$u11 - p$p4 * u$u12)
                     / sqrt(2))
        dp2db1 <- dn2 * (-p$p3 / sqrt(2)) * regr$X1
        dp2db2 <- dn2 * (-p$p4 / sqrt(2)) * regr$X2
        dp2db3 <- dn2 * (p$p5 / sqrt(2)) * regr$X3
        dp2db4 <- dn2 * (p$p6 / sqrt(2)) * regr$X4
        dp2dg2 <- dn2 * (u$u11 - u$u12) * dlink(u$u22 / sqrt(2)) * regr$Z2 / 2
        dp2dg4 <- dn2 * (u$u14 - u$u13) * dlink(u$u24 / sqrt(2)) * regr$Z4 / 2
        dp2 <- cbind(dp2db1, dp2db2, dp2db3, dp2db4, dp2dg2, dp2dg4)
        dp1 <- -dp2
    }

    dL1 <- (1 / p$p1) * dp1 + (1 / p$p3) * dp3
    dL2 <- (1 / p$p1) * dp1 + (1 / p$p4) * dp4
    dL3 <- (1 / p$p2) * dp2 + (1 / p$p5) * dp5
    dL4 <- (1 / p$p2) * dp2 + (1 / p$p6) * dp6

    ans <- matrix(NA, nrow = n, ncol = sum(rcols[1:6]))
    ans[y == 1, ] <- dL1[y == 1, ]
    ans[y == 2, ] <- dL2[y == 2, ]
    ans[y == 3, ] <- dL3[y == 3, ]
    ans[y == 4, ] <- dL4[y == 4, ]

    return(ans)
}

makeResponse122 <- function(yf)
{
    if (length(dim(yf))) {
        Y <- yf

        if (ncol(Y) > 3) {
            warning("only first three columns of response will be used")
            Y <- Y[, 1:3]
        }

        if (!all(unlist(Y) %in% c(0L, 1L)))
            stop("dummy responses must be dummy variables")

        if (ncol(Y) == 3) {
            ## this is for the case where y is specified as
            ##   (1's move) + (2's move if 1 moves L) + (2's move if 1 moves R)
            Y[, 2] <- ifelse(Y[, 1] == 1, Y[, 3], Y[, 2])
            ylevs <- c(paste("~", names(Y)[2], sep = ""),
                       names(Y)[2],
                       paste("~", names(Y)[3], sep = ""),
                       names(Y)[3])
        } else {
            ## this is for the case where y is specified as
            ##   (1's move) + (2's move)
            ylevs <- c(paste("~", names(Y)[1], ",~", names(Y)[2], sep = ""),
                       paste("~", names(Y)[1], ",", names(Y)[2], sep = ""),
                       paste(names(Y)[1], ",~", names(Y)[2], sep = ""),
                       paste(names(Y)[1], ",", names(Y)[2], sep = ""))
        }

        y <- numeric(nrow(Y))
        y[Y[, 1] == 0 & Y[, 2] == 0] <- 1
        y[Y[, 1] == 0 & Y[, 2] == 1] <- 2
        y[Y[, 1] == 1 & Y[, 2] == 0] <- 3
        y[Y[, 1] == 1 & Y[, 2] == 1] <- 4

        yf <- as.factor(y)
        levels(yf) <- ylevs
    } else {
        yf <- as.factor(yf)
        if (nlevels(yf) != 4) stop("dependent variable must have four values")
    }

    return(yf)
}

##' Strategic model with 2 players, 4 terminal nodes
##' 
##' Fits a strategic model with two players and four terminal nodes, as in the
##' game illustrated below in "Details".
##'
##' The model corresponds to the following extensive-form game:
##' \preformatted{
##' .        ___ 1 ___
##' .       /         \
##' .      /           \
##' .   2 /             \ 2
##' .    / \           / \
##' .   /   \         /   \
##' .  /     \       /     \
##' . u11    u12    u13    u14
##' . 0      u22    0      u24}
##'
##' For additional details on any of the function arguments or options, see
##' \code{\link{egame12}}.  The only difference is that the right-hand side of
##' \code{formulas} must have six components (rather than four) in this case.
##'
##' Ways to specify the dependent variable in \code{egame122}:
##' \itemize{
##' \item Numeric vector \code{y}, numbered 1 through 4, corresponding to the
##' outcomes as labeled in the game tree above.
##' \item Factor \code{y}, where \code{y} has four levels, corresponding in
##' order to the outcomes as labeled above.
##' \item Indicator variables \code{y1 + y2}, where \code{y1} indicates whether
##' Player 1 moves left or right, and \code{y2} indicates whether Player 2 moves
##' left or right.
##' \item Indicator variables \code{y1 + y2 + y3}, where \code{y1} indicates
##' whether Player 1 moves left or right, \code{y2} indicates Player 2's move in
##' case Player 1 moved left, and \code{y3} indicates Player 2's move in case
##' Player 1 moved right.  Non-observed values of \code{y2} and \code{y3} should
##' be set to \code{0}, \strong{not} \code{NA}, to ensure that observations are
##' not dropped when \code{na.action = na.omit}.}
##' @param formulas a list of six formulas, or a \code{Formula} object with six
##' right-hand sides.  See "Details" and "Examples".
##' @param data a data frame.
##' @param subset an optional logical vector specifying which observations from
##' \code{data} to use in fitting.
##' @param na.action how to deal with \code{NA}s in \code{data}.  Defaults to
##' the \code{na.action} setting of \code{\link{options}}.  See
##' \code{\link{na.omit}}
##' @param link whether to use a probit (default) or logit link structure,
##' @param type whether to use an agent-error ("agent", default) or
##' private-information ("private") stochastic structure.
##' @param startvals whether to calculate starting values for the optimization
##' from statistical backwards induction ("sbi", default), draw them from a
##' uniform distribution ("unif"), or to set them all to 0 ("zero")
##' @param fixedUtils numeric vector of values to fix for u11, u12, u13, u14,
##' u22, and u24.  \code{NULL} (the default) indicates that these should be
##' estimated with regressors, not fixed.
##' @param sdformula an optional list of formulas or a \code{\link{Formula}}
##' containing a regression equation for the scale parameter.  See
##' \code{\link{egame12}} for details.
##' @param sdByPlayer logical: if scale parameters are being estimated (i.e.,
##' \code{sdformula} or \code{fixedUtils} is non-\code{NULL}), should a separate
##' one be estimated for each player?  This option is ignored unless
##' \code{fixedUtils} or \code{sdformula} is specified.
##' @param boot integer: number of bootstrap iterations to perform (if any).
##' @param bootreport logical: whether to print status bar during bootstrapping.
##' @param profile output from running \code{\link{profile.game}} on a previous
##' fit of the model, used to generate starting values for refitting when an
##' earlier fit converged to a non-global maximum.
##' @param method character string specifying which optimization routine to use
##' (see \code{\link{maxLik}})
##' @param ... other arguments to pass to the fitting function (see
##' \code{\link{maxLik}}).
##' @return An object of class \code{c("game", "egame122")}.  See
##' \code{\link{egame12}} for a description of the \code{game} class.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com}) and Curtis
##' S. Signorino
##' @example inst/examples/egame122.r
egame122 <- function(formulas, data, subset, na.action,
                     link = c("probit", "logit"),
                     type = c("agent", "private"),
                     startvals = c("sbi", "unif", "zero"),
                     fixedUtils = NULL,
                     sdformula = NULL,
                     sdByPlayer = FALSE,
                     boot = 0,
                     bootreport = TRUE,
                     profile,
                     method = "BFGS",
                     ...)
{
    cl <- match.call()

    link <- match.arg(link)
    type <- match.arg(type)
    startvals <- match.arg(startvals)

    formulas <- checkFormulas(formulas)

    ## sanity checks
    if (!is.null(fixedUtils)) {
        if (length(fixedUtils) < 6)
            stop("fixedUtils must have 6 elements (u11, u12, u13, u14, u22, u24)")
        if (length(fixedUtils) > 6) {
            warning("only the first 6 elements of fixedUtils will be used")
            fixedUtils <- fixedUtils[1:6]
        }

        formulas <- update(formulas, . ~ 1 | 1 | 1 | 1 | 1 | 1)

        if (startvals == "sbi")
            startvals <- "zero"

        if (is.null(sdformula))
            sdformula <- if (sdByPlayer) Formula(~ 1 | 1) else Formula(~ 1)
    }

    if (!is.null(sdformula)) {
        sdformula <- checkFormulas(sdformula, argname = "sdformula")
        if (sdByPlayer && length(sdformula)[2] != 2)
            stop("`sdformula` should have two components (one for each player) on the right-hand side when sdByPlayer == TRUE")
        if (!sdByPlayer && length(sdformula)[2] != 1)
            stop("`sdformula` should have exactly one component on the right-hand side")
        formulas <- as.Formula(formula(formulas), formula(sdformula))
    }

    if (sdByPlayer && is.null(sdformula)) {
        warning("to estimate SDs by player, you must specify `sdformula` or `fixedUtils`")
        sdByPlayer <- FALSE
    }

    if (link == "logit" && type == "private") {
        warning("logit link cannot be used with private information model; changing to probit link")
        link <- "probit"
    }

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## make response variables (yf factor, y numeric)
    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    yf <- makeResponse122(yf)
    y <- as.numeric(yf)

    ## make list of regressor matrices
    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    rcols <- sapply(regr, ncol)

    ## starting values
    if (missing(profile) || is.null(profile)) {
        if (startvals == "zero") {
            sval <- rep(0, sum(rcols))
        } else if (startvals == "unif") {
            if (!hasArg(unif))
                unif <- c(-1, 1)
            sval <- runif(sum(rcols), unif[1], unif[2])
        } else {
            sval <- sbi122(y, regr, link)
            sval <- c(sval, rep(0, sum(rcols) - length(sval)))
        }
    } else {
        sval <- svalsFromProfile(profile)
    }

    ## identification check
    varNames <- lapply(regr, colnames)
    idCheck <- do.call(intersectAll, varNames[1:4])
    if (is.null(fixedUtils) && length(idCheck) > 0) {
        stop("Identification problem: the following variables appear in all four of player 1's utility equations: ",
             paste(idCheck, collapse = ", "))
    }

    ## make variable names
    prefixes <- paste(c(rep("u1(", 4), rep("u2(", 2)),
                      c(levels(yf), levels(yf)[2], levels(yf)[4]), ")",
                      sep = "")
    sdterms <- if (!is.null(sdformula)) { if (sdByPlayer) 2L else 1L } else 0L
    utils <- if (is.null(fixedUtils)) 1:6 else numeric(0)
    varNames <- makeVarNames(varNames, prefixes, utils, link, sdterms)
    hasColon <- varNames$hasColon
    names(sval) <- varNames$varNames

    ## use gradient iff no scale parameters being estimated
    gr <- if (is.null(sdformula)) logLikGrad122 else NULL

    ## deal with fixed utilities
    fvec <- rep(FALSE, length(sval))
    names(fvec) <- names(sval)
    if (!is.null(fixedUtils)) {
        sval[1:6] <- fixedUtils
        fvec[1:6] <- TRUE
    }

    results <- maxLik(logLik = logLik122, grad = gr, start = sval, fixed = fvec,
                      method = method, y = y, regr = regr, link = link, type =
                      type, ...)
    cc <- convergenceCriterion(method)
    if (!(results$code %in% cc)) {
        warning("Model fitting did not converge\nCode:", results$code,
                "\nMessage: ", results$message)
    }

    ## check local identification
    lid <- checkLocalID(results$hessian, fvec)
    if (!lid)
        warning("Hessian is not negative definite; coefficients may not be locally identified")

    if (boot > 0) {
        bootMatrix <-
            gameBoot(boot, report = bootreport, estimate = results$estimate, y =
                     y, regr = regr, fn = logLik122, gr = gr, fixed = fvec,
                     method = method, link = link, type = type, ...)
    }

    ## store output
    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getGameVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLik122(results$estimate, y = y, regr = regr, link = link, type = type)
    ans$call <- cl
    ans$convergence <- list(method = method, iter = nIter(results), code =
                            results$code, message = results$message, gradient =
                            !is.null(gr))
    ans$formulas <- formulas
    ans$link <- link
    ans$type <- type
    ans$model <- mf
    ans$xlevels <- .getXlevels(attr(mf, "terms"), mf)
    ans$y <- yf
    ans$equations <- names(hasColon)
    attr(ans$equations, "hasColon") <- hasColon
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    ans$localID <- lid

    class(ans) <- c("game", "egame122")

    return(ans)
}
