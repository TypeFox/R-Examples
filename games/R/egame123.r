##' @export
predict.egame123 <- function(object, newdata, type = c("outcome", "action"),
                             na.action = na.pass, ...)
{
    type <- match.arg(type)

    if (missing(newdata) || is.null(newdata)) {
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
    ans <- makeProbs123(object$coefficients, regr = regr, link = object$link,
                        type = object$type)

    if (type == "outcome") {
        ans <- data.frame(actionsToOutcomes123(ans, log.p = FALSE))
        names(ans) <- paste("Pr(", levels(object$y), ")", sep="")
    } else {
        ans <- as.data.frame(ans)
        names(ans)[1:2] <- paste("Pr(", c("", "~"), levels(object$y)[1], ")",
                                 sep="")
        names(ans)[3:4] <- paste("Pr(", c("", "~"), levels(object$y)[2], "|~",
                                 levels(object$y)[1], ")", sep="")
        names(ans)[5:6] <- paste("Pr(", levels(object$y)[3:4], "|~",
                                 levels(object$y)[1], ",~[", levels(object$y)[2],
                                 "])", sep = "")
        names(ans) <- gsub("~~", "", names(ans))
    }

    return(ans)
}

sbi123 <- function(y, regr, link)
{
    names(regr) <- character(length(regr))
    names(regr)[1:8] <- c("X1", "X3", "X5", "X6", "Z3", "Z5", "Z6", "W6")

    if (link == "probit") {
        fam <- binomial(link = "probit")
        linkfcn <- pnorm
    } else {
        fam <- binomial(link = "logit")
        linkfcn <- plogis
    }

    ## regression for player 3's choice
    reg3 <- regr$W6[y == 3 | y == 4, , drop = FALSE]
    y3 <- as.numeric(y == 4)[y == 3 | y == 4]
    m3 <- suppressWarnings(glm.fit(reg3, y3, family = fam))
    p6 <- as.numeric(regr$W6 %*% coef(m3))
    p6 <- linkfcn(p6)

    ## regression for player 2's choice
    reg2 <- cbind(-regr$Z3, (1-p6) * regr$Z5, p6 * regr$Z6)
    reg22 <- reg2[y != 1, , drop = FALSE]
    y22 <- as.numeric(y != 2)[y != 1]
    m22 <- suppressWarnings(glm.fit(reg22, y22, family = fam))
    p4 <- as.numeric(reg2 %*% coef(m22))
    p4 <- linkfcn(p4)

    ## regression for player 1's choice
    reg1 <- cbind(-regr$X1, (1-p4) * regr$X3, p4 * (1-p6) * regr$X5,
                  p4 * p6 * regr$X6)
    y1 <- as.numeric(y != 1)
    m1 <- suppressWarnings(glm.fit(reg1, y1, family = fam))

    ans <- sqrt(2) * c(coef(m1), coef(m22), coef(m3))
    return(ans)
}

makeSDs123 <- function(b, regr, type)
{
    sds <- vector("list", if (type == "private") 9L else 6L)
    regr <- regr[-(1:8)]
    rcols <- sapply(regr, ncol)

    if (length(rcols) == 1L) {  ## sdByPlayer == FALSE
        v <- exp(as.numeric(regr[[1]] %*% b))
        for (i in 1:length(sds)) sds[[i]] <- v
    } else {
        b1 <- b[1:rcols[1]]
        b2 <- b[(rcols[1]+1):(rcols[1]+rcols[2])]
        b3 <- b[(rcols[1]+rcols[2]+1):length(b)]
        v1 <- exp(as.numeric(regr[[1]] %*% b1))
        v2 <- exp(as.numeric(regr[[2]] %*% b2))
        v3 <- exp(as.numeric(regr[[3]] %*% b3))

        if (type == "private") {
            sds[[1]] <- sds[[2]] <- sds[[3]] <- sds[[4]] <- v1
            sds[[5]] <- sds[[6]] <- sds[[7]] <- v2
            sds[[8]] <- sds[[9]] <- v3
        } else {
            sds[[1]] <- sds[[2]] <- v1
            sds[[3]] <- sds[[4]] <- v2
            sds[[5]] <- sds[[6]] <- v3
        }
    }

    return(sds)
}

makeProbs123 <- function(b, regr, link, type)
{
    utils <- makeUtils(b, regr, nutils = 8,
                       unames = c("u11", "u13", "u15", "u16", "u23", "u25",
                       "u26", "u36"))

    if (length(utils$b) == 0) {  ## variance unparameterized
        sds <- as.list(rep(1, 9))
    } else {
        sds <- makeSDs123(utils$b, regr, type)
    }

    linkfcn <- switch(link,
                      logit = function(x, sd = 1) plogis(x, scale = sd),
                      probit = pnorm)

    if (type == "private") {
        sd6 <- sqrt(sds[[8]]^2 + sds[[9]]^2)
    } else {
        sd6 <- sqrt(sds[[5]]^2 + sds[[6]]^2)
    }
    p6 <- finiteProbs(linkfcn(utils$u36, sd = sd6))
    p5 <- 1 - p6

    if (type == "private") {
        sd4 <- sqrt(p5^2 * sds[[6]]^2 + p6^2 * sds[[7]]^2 + sds[[5]]^2)
    } else {
        sd4 <- sqrt(sds[[3]]^2 + sds[[4]]^2)
    }
    p4 <- p5 * utils$u25 + p6 * utils$u26 - utils$u23
    p4 <- finiteProbs(linkfcn(p4, sd = sd4))
    p3 <- 1 - p4

    if (type == "private") {
        sd2 <- sqrt(p3^2 * sds[[2]]^2 + p4^2 * p5^2 * sds[[3]]^2 +
                    p4^2 * p6^2 * sds[[4]]^2 + sds[[1]]^2)
    } else {
        sd2 <- sqrt(sds[[1]]^2 + sds[[2]]^2)
    }
    p2 <- p3 * utils$u13 + p4 * (p5 * utils$u15 + p6 * utils$u16) - utils$u11
    p2 <- finiteProbs(linkfcn(p2, sd = sd2))
    p1 <- 1 - p2

    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}

actionsToOutcomes123 <- function(probs, log.p = TRUE)
{
    probs <- log(do.call(cbind, probs))
    ans <- cbind(probs[, 1],
                 probs[, 2] + probs[, 3],
                 probs[, 2] + probs[, 4] + probs[, 5],
                 probs[, 2] + probs[, 4] + probs[, 6])

    if (!log.p) ans <- exp(ans)
    return(ans)
}

logLik123 <- function(b, y, regr, link, type, ...)
{
    probs <- makeProbs123(b, regr, link, type)
    logProbs <- actionsToOutcomes123(probs, log.p = TRUE)
    ans <- logProbs[cbind(1:nrow(logProbs), y)]
    return(ans)
}

logLikGrad123 <- function(b, y, regr, link, type, ...)
{
    names(regr) <- character(length(regr))
    names(regr)[1:8] <- c("X1", "X3", "X5", "X6", "Z3", "Z5", "Z6", "W6")

    u <- makeUtils(b, regr, nutils = 8,
                   unames = c("u11", "u13", "u15", "u16", "u23", "u25",
                   "u26", "u36"))
    p <- makeProbs123(b, regr, link, type)
    eu24 <- p$p5 * u$u25 + p$p6 * u$u26 - u$u23
    eu12 <- p$p3 * u$u13 + p$p4*(p$p5*u$u15 + p$p6*u$u16) - u$u11
    eu14c2 <- p$p5*u$u15 + p$p6*u$u16 - u$u13
    rcols <- sapply(regr, ncol)
    n <- nrow(regr$X1)

    if (link == "probit" && type == "private") {
        dp6db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp6dg <- matrix(0L, nrow = n, ncol = sum(rcols[5:7]))
        dp6du <- dnorm(u$u36 / sqrt(2)) * regr$W6 / sqrt(2)
        dp6 <- cbind(dp6db, dp6dg, dp6du)
        dp5 <- -dp6

        denom4 <- sqrt(1 + p$p5^2 + p$p6^2)
        phi24 <- dnorm(eu24 / denom4)
        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg3 <- -phi24 * regr$Z3 / denom4
        dp4dg5 <- p$p5 * phi24 * regr$Z5 / denom4
        dp4dg6 <- p$p6 * phi24 * regr$Z6 / denom4
        dp4dg <- cbind(dp4dg3, dp4dg5, dp4dg6)
        dp4du <- phi24 * ((u$u26-u$u25)*denom4 - eu24*(p$p6-p$p5)/denom4)
        dp4du <- (dp4du / denom4^2) * dp6du
        dp4 <- cbind(dp4db, dp4dg, dp4du)
        dp3 <- -dp4

        denom2 <- sqrt(1 + p$p3^2 + (p$p4^2)*(p$p5^2) + (p$p4^2)*(p$p6^2))
        phi12 <- dnorm(eu12 / denom2)
        dp2db1 <- -phi12 * regr$X1 / denom2
        dp2db3 <- p$p3 * phi12 * regr$X3 / denom2
        dp2db5 <- p$p4 * p$p5 * phi12 * regr$X5 / denom2
        dp2db6 <- p$p4 * p$p6 * phi12 * regr$X6 / denom2
        dp2dg <- (eu14c2 * denom2 - eu12*(p$p4*(p$p5^2+p$p6^2)-p$p3)/denom2)
        dp2dg <- phi12 * (dp2dg / denom2^2) * dp4dg
        deu12du <- u$u13*(-dp4du) + u$u15*(p$p4*(-dp6du) + p$p5*dp4du) +
            u$u16*(p$p4*dp6du + p$p6*dp4du)
        ddenom2du <- p$p3*(-dp4du) + (p$p4^2)*p$p5*(-dp6du) +
            p$p4*(p$p5^2)*dp4du + (p$p4^2)*p$p6*dp6du + p$p4*(p$p6^2)*dp4du
        dp2du <- deu12du * denom2 - eu12 * ddenom2du / denom2
        dp2du <- phi12 * (dp2du / denom2^2)
        dp2 <- cbind(dp2db1, dp2db3, dp2db5, dp2db6, dp2dg, dp2du)
        dp1 <- -dp2
    } else if (type == "agent") {
        dlink <- switch(link,
                        logit = dlogis,
                        probit = dnorm)

        phi12 <- dlink(eu12 / sqrt(2))
        phi24 <- dlink(eu24 / sqrt(2))
        phi36 <- dlink(u$u36 / sqrt(2))

        dp6db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp6dg <- matrix(0L, nrow = n, ncol = sum(rcols[5:7]))
        dp6du <- phi36 * regr$W6 / sqrt(2)
        dp6 <- cbind(dp6db, dp6dg, dp6du)
        dp5 <- -dp6

        dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:4]))
        dp4dg3 <- -phi24 * regr$Z3 / sqrt(2)
        dp4dg5 <- p$p5 * phi24 * regr$Z5 / sqrt(2)
        dp4dg6 <- p$p6 * phi24 * regr$Z6 / sqrt(2)
        dp4du <- phi24 * phi36 * (u$u26 - u$u25) * regr$W6 / 2
        dp4 <- cbind(dp4db, dp4dg3, dp4dg5, dp4dg6, dp4du)
        dp3 <- -dp4

        dp2db1 <- -phi12 * regr$X1 / sqrt(2)
        dp2db3 <- p$p3 * phi12 * regr$X3 / sqrt(2)
        dp2db5 <- p$p4 * p$p5 * phi12 * regr$X5 / sqrt(2)
        dp2db6 <- p$p4 * p$p6 * phi12 * regr$X6 / sqrt(2)
        dp2dg3 <- -eu14c2 * phi12 * phi24 * regr$Z3 / 2
        dp2dg5 <- p$p5 * eu14c2 * phi12 * phi24 * regr$Z5 / 2
        dp2dg6 <- p$p6 * eu14c2 * phi12 * phi24 * regr$Z6 / 2
        dp2du <- -phi12 * phi36 *
            (p$p4*(u$u15-u$u16)*sqrt(2) + eu14c2*(u$u25-u$u26)*phi24) *
                regr$W6 / (2*sqrt(2))
        dp2 <- cbind(dp2db1, dp2db3, dp2db5, dp2db6, dp2dg3, dp2dg5, dp2dg6,
                     dp2du)
        dp1 <- -dp2
    }

    dL1 <- (1 / p$p1) * dp1
    dL2 <- (1 / p$p2) * dp2 + (1 / p$p3) * dp3
    dL3 <- (1 / p$p2) * dp2 + (1 / p$p4) * dp4 + (1 / p$p5) * dp5
    dL4 <- (1 / p$p2) * dp2 + (1 / p$p4) * dp4 + (1 / p$p6) * dp6

    ans <- matrix(NA, nrow = n, ncol = sum(rcols[1:8]))
    ans[y == 1, ] <- dL1[y == 1, ]
    ans[y == 2, ] <- dL2[y == 2, ]
    ans[y == 3, ] <- dL3[y == 3, ]
    ans[y == 4, ] <- dL4[y == 4, ]

    return(ans)
}

makeResponse123 <- function(yf)
{
    if (length(dim(yf))) {  ## yf is a matrix of dummies
        if (ncol(yf) == 2) {
            stop("response must be specified as a single vector or three dummy variables")
        } else if (ncol(yf) > 3) {
            warning("only first three columns of response will be used")
            yf <- yf[, 1:3]
        }

        if (!(all(unlist(yf) %in% c(0L, 1L))))
            stop("dummy responses must be dummy variables")

        ylevs <- c(paste("~", names(yf)[1], sep = ""),
                   paste(names(yf)[1], ",~", names(yf)[2], sep = ""),
                   paste(names(yf)[1], ",", names(yf)[2], ",~", names(yf)[3],
                         sep = ""),
                   paste(names(yf)[1], names(yf)[2], names(yf)[3], sep = ","))

        y <- integer(nrow(yf))
        y[yf[, 1] == 0] <- 1L
        y[yf[, 1] == 1 & yf[, 2] == 0] <- 2L
        y[yf[, 1] == 1 & yf[, 2] == 1 & yf[, 3] == 0] <- 3L
        y[yf[, 1] == 1 & yf[, 2] == 1 & yf[, 3] == 1] <- 4L
        yf <- as.factor(y)
        levels(yf) <- ylevs
    } else {                ## yf is a vector
        yf <- as.factor(yf)
        if (nlevels(yf) != 4) stop("dependent variable must have four values")
    }

    return(yf)
}

##' Strategic model with 3 players, 4 terminal nodes
##' 
##' Fits a strategic model with three players and four terminal nodes, as in the
##' game illustrated below in "Details".
##'
##' The model corresponds to the following extensive-form game:
##' \preformatted{
##' .     1
##' .     /\
##' .    /  \
##' .   /    \ 2
##' .  u11   /\
##' .       /  \
##' .      /    \
##' .    u13     \ 3
##' .    u23     /\
##' .           /  \
##' .          /    \
##' .         u15   u16
##' .         u25   u26
##' .         0     u36}
##'
##' For additional details on any of the function arguments or options, see
##' \code{\link{egame12}}.  The only difference is that the right-hand side of
##' \code{formulas} must have eight components (rather than four) in this case.
##'
##' Ways to specify the dependent variable in \code{egame123}:
##' \itemize{
##' \item Numeric vector \code{y} containing 4 unique values, corresponding to
##' the outcomes (in order from left to right) as labeled in the game tree
##' above.
##' \item Factor \code{y}, where \code{y} has four levels, corresponding in
##' order to the outcomes as labeled above.
##' \item Indicator variables \code{y1 + y2 + y3}, where \code{y1} indicates
##' whether Player 1 moves left or right, \code{y2} indicates Player 2's move,
##' and \code{y3} indicates Player 3's move.  Non-observed values of \code{y2}
##' and \code{y3} (where the game ended before the move could be made) should be
##' set to \code{0}, \strong{not} \code{NA}, to ensure that observations are not
##' dropped when \code{na.action = na.omit}.}
##' @param formulas a list of eight formulas, or a \code{\link{Formula}} object
##' with eight right-hand sides.  See "Details" and "Examples".
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
##' @param fixedUtils numeric vector of values to fix for u11, u13, u15, u16,
##' u23, u25, u26, and u36.  \code{NULL} (the default) indicates that these
##' should be estimated with regressors, not fixed.
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
##'ode{\link{maxLik}}).
##' @return An object of class \code{c("game", "egame123")}.  See
##' \code{\link{egame12}} for a description of the \code{game} class.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @example inst/examples/egame123.r
egame123 <- function(formulas, data, subset, na.action,
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

    if (!is.null(fixedUtils)) {  ## error checking for fixed utilities
        if (length(fixedUtils) < 8)
            stop("fixedUtils must have 8 elements (u11, u13, u15, u16, u23, u25, u26, u36)")
        if (length(fixedUtils) > 8) {
            warning("only the first 8 elements of fixedUtils will be used")
            fixedUtils <- fixedUtils[1:8]
        }

        formulas <- update(formulas, . ~ 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1)

        if (startvals == "sbi")
            startvals <- "zero"

        if (is.null(sdformula))
            sdformula <- if (sdByPlayer) Formula(~ 1 | 1 | 1) else Formula(~ 1)
    }

    if (!is.null(sdformula)) {  ## error checking for parameterized variance
        sdformula <- checkFormulas(sdformula, argname = "sdformula")
        if (sdByPlayer && length(sdformula)[2] != 3)
            stop("'sdformula' should have three components (one for each player) on the right-hand side when sdByPlayer == TRUE")
        if (!sdByPlayer && length(sdformula)[2] != 1)
            stop("'sdformula' should have exactly one component on the right-hand side")

        ## make one big Formula object with all utility and variance equations
        ## on the right-hand side
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

    yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    yf <- makeResponse123(yf)
    y <- as.numeric(yf)

    regr <- list()
    for (i in seq_len(length(formulas)[2]))
        regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
    rcols <- sapply(regr, ncol)

    ## calculate starting values
    if (missing(profile) || is.null(profile)) {
        if (startvals == "zero") {
            sval <- rep(0, sum(rcols))
        } else if (startvals == "unif") {
            if (!hasArg(unif))
                unif <- c(-1, 1)
            sval <- runif(sum(rcols), unif[1], unif[2])
        } else {
            sval <- sbi123(y, regr, link)
            sval <- c(sval, rep(0, sum(rcols) - length(sval)))
        }
    } else {
        sval <- svalsFromProfile(profile)
    }

    ## identification check
    varNames <- lapply(regr, colnames)
    idCheck <- do.call(intersectAll, varNames[1:4])
    idCheck2 <- do.call(intersectAll, varNames[5:7])
    if (is.null(fixedUtils) && (length(idCheck) > 0)) {
        stop("Identification problem: the following variables appear in all four of player 1's utility equations: ",
             paste(idCheck, collapse =", "))
    } else if (is.null(fixedUtils) && length(idCheck2 > 0)) {
        stop("Identification problem: the following variables appear in all three of player 2's utility equations: ",
             paste(idCheck2, collapse =", "))
    }

    ## variable naming
    prefixes <- paste(c(rep("u1(", 4), rep("u2(", 3), "u3("),
                      c(levels(yf), levels(yf)[2:4], levels(yf)[4]), ")",
                      sep = "")
    sdterms <- if (!is.null(sdformula)) { if (sdByPlayer) 3L else 1L } else 0L
    utils <- if (is.null(fixedUtils)) 1:8 else numeric(0)
    varNames <- makeVarNames(varNames, prefixes, utils, link, sdterms)
    hasColon <- varNames$hasColon
    names(sval) <- varNames$varNames

    ## use gradient only if variance isn't parameterized
    gr <- if (is.null(sdformula)) logLikGrad123 else NULL

    fvec <- rep(FALSE, length(sval))
    names(fvec) <- names(sval)
    if (!is.null(fixedUtils)) {
        sval[1:8] <- fixedUtils
        fvec[1:8] <- TRUE
    }

    results <- maxLik(logLik = logLik123, grad = gr, start = sval, fixed = fvec,
                      method = method, y = y, regr = regr, link = link, type =
                      type, ...)

    ## check for convergence
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
                     y, regr = regr, fn = logLik123, gr = gr, fixed = fvec,
                     method = method, link = link, type = type, ...)
    }

    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getGameVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLik123(results$estimate, y = y, regr = regr, link = link, type =
                  type)
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
    ans$equations <- structure(names(hasColon), hasColon = hasColon)
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    ans$localID <- lid

    class(ans) <- c("game", "egame123")

    return(ans)
}
