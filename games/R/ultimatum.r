##' Lambert's W
##' 
##' Solves for W in the equation \eqn{W e^W = x}{W * exp(W) = x}.
##'
##' The function is based on the code given in Barry et al. (1995).  It is used
##' to calculate fitted values for the \code{\link{ultimatum}} model.
##'
##' If negative values of \code{x} are supplied, \code{NaN}s will likely be
##' returned.
##' @param x vector of values to solve for.
##' @return Solutions to Lambert's W for each value in \code{x}.
##' @export
##' @references D.A. Barry, P.J. Culligan-Hensley, and S.J. Barry.  1995.  "Real
##' Values of the W-Function."  \emph{ACM Transactions on Mathematical Software}
##' 21(2):161--171.
##' @author Curt Signorino (\email{curt.signorino@@rochester.edu})
##' @example inst/examples/LW.r
LW <- function(x)
{
    ## Note: there is a Lambert's W package, but it depends on the gsl (GNU
    ## Scientific Library) package, which is hellish for Windows users to
    ## install, hence our hand-rolled Lambert's W.
    
    eW <- function(x, W)
    {
        zn <- log(x/W)-W
        first <- zn/(1+W)
        common <- 2*(1+W)*(1+W+2*zn/3)
        first*(common-zn)/(common-2*zn)
    }
    
    lx <- log(x)
    wlt <- (x*(1+4/3*x))/(1+x*(7/3+5/6*x))
    wgt <- lx-(24*(lx*(lx+2)-3))/(lx*(7*lx+58)+127)
    W1 <- ifelse(x < 0.7385, wlt, wgt)
    
    xge20 <- x >= 20
    i <- seq_along(x)
    ixge20 <- i[xge20]
    a1 <- 1.124491989777808
    b1 <- .4225028202459761
    xg <- x[ixge20]
    h <- exp(-a1/(b1+log(xg)))
    Wp2 <- log(xg/log(xg/(log(xg))^h))  # for x>20
    W1[ixge20] <- Wp2;                   

    ## iteration for improved accuracy
    W2 <- W1*(1+eW(x,W1))
    W3 <- W2*(1+eW(x,W2))
    W4 <- W3*(1+eW(x,W3))               # enough

    return(W4)
}

##' @export
predict.ultimatum <- function(object, newdata, na.action = na.pass,
                              n.sim = 1000, ...)
{
    if (missing(newdata)) {
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

    X <- model.matrix(object$formulas, data = mf, rhs = 1)
    Z <- model.matrix(object$formulas, data = mf, rhs = 2)

    ## extract relevant model components
    b <- object$coefficients
    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(Z))
    b <- head(b, length(b) - ncol(Z))
    fit1 <- as.numeric(X %*% b)
    fit2 <- as.numeric(Z %*% g)
    maxOffer <- object$maxOffer
    minOffer <- object$minOffer

    ## There isn't an analytic expression for the expected value of Player 1's
    ## optimal offer, so we need to simulate
    e1 <- matrix(rlogis(nrow(X) * n.sim, scale = s1), nrow = nrow(X))
    bestOffer <- (maxOffer - fit1 - e1 - s2 - fit2) / s2
    bestOffer <- exp(bestOffer)
    bestOffer[bestOffer > .Machine$double.xmax] <- .Machine$double.xmax
    bestOffer[bestOffer < .Machine$double.xmin] <- .Machine$double.xmin
    bestOffer <- maxOffer - fit1 - e1 - s2 * (1 + LW(bestOffer))
    offer <- bestOffer
    offer[offer < minOffer] <- minOffer
    offer[offer > maxOffer] <- maxOffer

    e2 <- matrix(rlogis(nrow(X) * n.sim, scale = s2), nrow = nrow(X))
    accept <- offer > fit2 + e2

    ## expected offer and probability of acceptance
    Ey <- rowMeans(offer)
    PrA <- rowMeans(accept)

    ans <- as.data.frame(cbind(Ey, PrA))
    names(ans) <- c("E(offer)", "Pr(accept)")
    return(ans)
}

##
## INPUT:
## x: numeric vector
##
## RETURN:
## numeric vector, replacing Inf with largest representable values
##
finitize <- function(x)
{
    x <- ifelse(is.finite(x), x, sign(x) * .Machine$double.xmax)
    return(x)
}

offerCDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    1 / (1 + exp((maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1))
}

offerPDF <- function(y, maxOffer, fit1, fit2, s1, s2)
{
    num <- exp(-(maxOffer - y - fit1 - s2 * (1 + exp((y - fit2) / s2))) / s1)
    denom <- s1 * (1 + num)^2
    ans <- (num / denom) * (1 + exp((y - fit2) / s2))
    return(ans)
}

logLikUlt <- function(b, y, acc, regr, minOffer=0, maxOffer, offerOnly, offertol, ...)
{
    ## extract parameters
    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(regr$Z))
    b <- head(b, length(b) - ncol(regr$Z))

    ## fitted reservation values
    fit1 <- as.numeric(regr$X %*% b)
    fit2 <- as.numeric(regr$Z %*% g)

    ## probability of acceptance (using finiteProbs to avoid numerical issues)
    prAccept <- finiteProbs(1 / (1 + exp(-(y - fit2) / s2)))

    ## probability of making an offer of 0 ('lowball') or the maximal offer
    ## ('highball'), and interior density of observed offer
    lowball <- finiteProbs(offerCDF(minOffer, maxOffer, fit1, fit2, s1, s2))
    highball <- finiteProbs(1 - offerCDF(maxOffer, maxOffer, fit1, fit2, s1,
                                         s2))
    interior <- offerPDF(y, maxOffer, fit1, fit2, s1, s2)

    ## identify maximal/minimal offers
    isMax <- abs(y - maxOffer) < offertol
    isMin <- abs(y - minOffer) < offertol

    ## evaulate log-likelihood for offers
    ans1 <- ifelse(isMax, highball, ifelse(isMin, lowball, interior))
    ans1 <- replace(ans1, ans1 < .Machine$double.eps, .Machine$double.eps)
    ans <- ans1
    attr(ans, "offer") <- log(ans1)

    ## evaluate log-likelihood for acceptance (if observed)
    if (!offerOnly) {
        ans2 <- finiteProbs(ifelse(acc == 1, prAccept, 1 - prAccept))
        ans <- ans * ans2
        attr(ans, "accept") <- log(ans2)
    }

    ans <- log(ans)
    return(ans)
}

logLikGradUlt <- function(b, y, acc, regr, minOffer=0, maxOffer, offerOnly, offertol, ...)
{
    isMax <- abs(y - maxOffer) < offertol
    isMin <- abs(y - minOffer) < offertol

    s1 <- exp(b[length(b) - 1])
    s2 <- exp(b[length(b)])
    b <- head(b, length(b) - 2)
    g <- tail(b, ncol(regr$Z))
    b <- head(b, length(b) - ncol(regr$Z))

    ## 'finitize' used liberally here to avoid numerical problems -- experience
    ## suggests it would be ill-advised to remove it, since all the exponentials
    ## can blow up, and ratios of them result in NaNs, which crash the fitting
    ## procedure
    fit1 <- finitize(as.numeric(regr$X %*% b))
    fit2 <- finitize(as.numeric(regr$Z %*% g))

    ey <- finitize(exp((y - fit2) / s2))
    ey1ey <- ey / (1 + ey)
    Qy <- finitize((-maxOffer + y + fit1 + s2 * (1 + ey)) / s1)
    Qy1Qy <- finitize(exp(Qy)) / finitize(1 + exp(Qy))
    Qy2Qy <- finitize(exp(-Qy)) / finitize(1 + exp(-Qy))

    dfdb <- ((1 - 2 * Qy1Qy) / s1) * regr$X
    dFdb <- Qy2Qy * (regr$X / s1)
    d1mFdb <- dFdb - regr$X / s1

    dfdg <- ((-ey1ey / s2) + (-1 + 2 * Qy1Qy) * (ey / s1)) * regr$Z
    dFdg <- Qy2Qy * ey * (-regr$Z / s1)
    d1mFdg <- dFdg + (ey / s1) * regr$Z

    dfdlns1 <- 2 * Qy1Qy * Qy - Qy - 1
    dFdlns1 <- Qy2Qy * (-Qy)
    d1mFdlns1 <- dFdlns1 + Qy

    dfdlns2 <- ((fit2 - y) * ey + s2 * (1 + ey)) / s1
    dfdlns2 <- dfdlns2 - 2 * Qy1Qy * dfdlns2 + ey1ey * ((fit2 - y) / s2)
    dFdlns2 <- (Qy2Qy / s1) * (ey * (fit2 - y) + s2 * (1 + ey))
    d1mFdlns2 <- dFdlns2 - (ey * (fit2 - y) + s2 * (1 + ey)) / s1

    df <- cbind(dfdb, dfdg, dfdlns1, dfdlns2)
    dF <- cbind(dFdb, dFdg, dFdlns1, dFdlns2)
    d1mF <- cbind(d1mFdb, d1mFdg, d1mFdlns1, d1mFdlns2)

    ans <- df
    ans[isMin, ] <- dF[isMin, ]
    ans[isMax, ] <- d1mF[isMax, ]

    ## gradient contribution in case acceptance is observed
    if (!offerOnly) {
        ey2ey <-
            finitize(exp((fit2 - y) / s2)) / finitize(1 + exp((fit2 - y) / s2))

        dPdb <- matrix(0L, nrow = nrow(regr$X), ncol = ncol(regr$X))
        dPdg <- -ey2ey * (regr$Z / s2)
        dPdlns1 <- 0L
        dPdlns2 <- -ey2ey * ((y - fit2) / s2)

        d1mPdb <- dPdb
        d1mPdg <- regr$Z / s2 + dPdg
        d1mPdlns1 <- 0L
        d1mPdlns2 <- ((y - fit2) / s2) + dPdlns2

        dAcc <- cbind(dPdb, dPdg, dPdlns1, dPdlns2)
        dRej <- cbind(d1mPdb, d1mPdg, d1mPdlns1, d1mPdlns2)

        ans2 <- dAcc
        ans2[acc == 0, ] <- dRej[acc == 0, ]

        ans <- ans + ans2
    }

    return(ans)
}

##' Statistical ultimatum game
##' 
##' Estimates the statistical ultimatum game described in Ramsay and Signorino
##' (2009), illustrated below in "Details".
##'
##' The model corresponds to the following extensive-form game, described in
##' Ramsay and Signorino (2009):
##' \preformatted{
##' .       1
##' .      / \
##' .     /   \
##' .    /     \ y in [0, Q]
##' .   /       \
##' .   ---------
##' .       /\  2
##' .      /  \
##' .     /    \
##' .    /      \
##' . Q - y     R1
##' . y         R2}
##' Q refers to the maximum feasible offer (the argument \code{maxOffer}).
##'
##' The two equations on the right-hand side of \code{formulas} refer to Player
##' 1's and Player 2's reservation values respectively.  The left-hand side
##' should take the form \code{offer + acceptance}, where \code{outcome}
##' contains the numeric value of the offer made and \code{acceptance} is an
##' indicator for whether it was accepted.  (If \code{outcome} is set to
##' "offer", the acceptance indicator can be omitted.  See below for more.)
##'
##' The \code{outcome} argument refers to whether the outcome of interest is
##' just the level of the offer made, or both the level of the offer and whether
##' it was accepted.  If acceptance was unobserved, then \code{outcome} should
##' be set to "offer".  If so, the estimates for Player 2's reservation value
##' should be interpreted as Player 1's expectations about these parameters.  It
##' may also be useful to set \code{outcome} to "offer" even if acceptance data
##' are available, for the purpose of comparing the strategic model to other
##' models of offer levels (as in Ramsay and Signorino 2009).  If an acceptance
##' variable is specified but \code{outcome} is set to "offer", the acceptance
##' data will be used for starting values but not in the actual fitting.
##'
##' Numerical instability is not uncommon in the statistical ultimatum game,
##' especially when the scale parameters are being estimated.
##' @param formulas a list of two formulas, or a \code{\link{Formula}} object
##' with two right-hand sides.  See "Details" and the examples below.
##' @param data data frame containing the variables in the model.
##' @param subset optional logical expression specifying which observations from
##' \code{data} to use in fitting.
##' @param na.action how to deal with \code{NA}s in \code{data}.  Defaults to
##' the \code{na.action} setting of \code{\link{options}}.  See
##' \code{\link{na.omit}}.
##' @param minOffer numeric: the lowest offer Player 1 could feasibly make
##' (default 0).
##' @param maxOffer numeric: the highest offer Player 1 could feasibly make.
##' @param offertol numeric: offers within \code{offertol} of
##' \code{minOffer}/\code{maxOffer} will be considered to be at the
##' minimum/maximum.  (This is used to prevent floating-point problems and
##' need not be changed in most applications.)
##' @param s1 numeric: scale parameter for Player 1.  If \code{NULL} (the
##' default), the parameter will be estimated.
##' @param s2 numeric: scale parameter for Player 2.  If \code{NULL} (the
##' default), the parameter will be estimated.
##' @param outcome the outcome of interest: just Player 1's offer ("offer") or
##' both the offer and its acceptance ("both").  See "Details".
##' @param boot integer: number of bootstrap iterations to perform (if any).
##' @param bootreport logical: whether to print status bar when performing
##' bootstrap iterations.
##' @param profile output from running \code{\link{profile.game}} on a previous
##' fit of the model, used to generate starting values for refitting when an
##' earlier fit converged to a non-global maximum.
##' @param method character string specifying which optimization routine to use
##' (see \code{\link{maxLik}})
##' @param ... other arguments to pass to the fitting function (see
##' \code{\link{maxLik}}).
##' @param reltol numeric: relative convergence tolerance level (see
##' \code{\link{optim}}).  Use of values higher than the default is discouraged.
##' @return An object of class \code{c("game", "ultimatum")}.  For details on
##' the \code{game} class, see \code{\link{egame12}}.  The \code{ultimatum}
##' class is just for use in the generation of predicted values (see
##' \code{\link{predProbs}}) and profiling (see \code{\link{profile.game}}).
##' @export
##' @references Kristopher W. Ramsay and Curtis S. Signorino.  2009.  "A
##' Statistical Model of the Ultimatum Game."  Available online at
##' \url{http://www.rochester.edu/college/psc/signorino/research/RamsaySignorino_Ultimatum.pdf}.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com}) and Curtis
##' S. Signorino
##' @example inst/examples/ultimatum.r
ultimatum <- function(formulas, data, subset, na.action, minOffer=0,
                      maxOffer, offertol = sqrt(.Machine$double.eps),
                      s1 = NULL, s2 = NULL,
                      outcome = c("both", "offer"),
                      boot = 0,
                      bootreport = TRUE,
                      profile,
                      method = "BFGS",
                      ...,
                      reltol = 1e-12)
{
    cl <- match.call()

    outcome <- match.arg(outcome)
    offerOnly <- outcome == "offer"
    s1null <- is.null(s1)
    s2null <- is.null(s2)

    formulas <- checkFormulas(formulas)

    ## make the model frame
    mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formulas
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## sanity checks/determine whether acceptance was specified
    ya <- model.part(formulas, mf, lhs = 1, drop = TRUE)
    if (length(dim(ya))) {
        y <- ya[, 1]
        a <- ya[, 2]
        if (!all(unique(a) %in% c(0, 1)))
            stop("acceptance variable must be binary")
    } else {
        if (!offerOnly) {
            stop("a dependent variable for acceptance must be specified if",
                 " `outcome == \"both\"`; see `?ultimatum`")
        }
        y <- ya
        a <- NULL
    }
    if (any(y > maxOffer))
        stop("observed offers greater than maxOffer")
	if (sum(y<minOffer)>0)
		stop("observed offer less than minOffer")
		
    ## retrieve the regressors for each player's reservation value
    regr <- list()
    regr$X <- model.matrix(formulas, data = mf, rhs = 1)
    regr$Z <- model.matrix(formulas, data = mf, rhs = 2)

    ## scale terms
    s1 <- if (s1null) log(sd(y)) else log(s1)
    s2 <- if (s2null) log(sd(y)) else log(s2)

    if (missing(profile) || is.null(profile)) {
        ## run a kind of pseudo-sbi to get starting values (performance of this
        ## is iffy, especially for getting the probability of acceptance)
        aa <- if (!is.null(a)) a else as.numeric(y >= mean(y))
        m2 <- suppressWarnings(glm.fit(regr$Z, aa,
                                       family = binomial(link = "logit"),
                                       intercept = FALSE, offset =
                                       as.numeric(-y)))
        m1 <- lsfit(regr$X, maxOffer - y, intercept = FALSE)
        sval <- c(m1$coefficients, m2$coefficients, s1, s2)

        ## see if the pseudo-sbi yields a finite log-likelihood; if not, just
        ## set everything to 0 except the constants
        firstTry <- logLikUlt(sval, y = y, acc = a, regr = regr, minOffer=minOffer,
        					  maxOffer = maxOffer, offerOnly = offerOnly, offertol =
                              offertol)
        if (!is.finite(sum(firstTry))) {
            sval <- c(maxOffer - mean(y), rep(0, length(m1$coefficients) - 1),
                      maxOffer - mean(y), rep(0, length(m2$coefficients) - 1),
                      s1, s2)
        }
    } else {
        sval <- svalsFromProfile(profile)
    }

    ## give names to starting values
    names(sval) <- c(paste("R1", colnames(regr$X), sep = ":"),
                     paste("R2", colnames(regr$Z), sep = ":"),
                     "log(s1)", "log(s2)")

    ## determine which parameters are to be fixed
    fvec <- logical(length(sval))
    if (!s1null) fvec[length(fvec) - 1] <- TRUE
    if (!s2null) fvec[length(fvec)] <- TRUE
    names(fvec) <- names(sval)

    results <- maxLik(logLik = logLikUlt, grad = logLikGradUlt, start = sval,
                      fixed = fvec, method = method, y = y, acc = a, regr =
                      regr, minOffer = minOffer, maxOffer = maxOffer, 
                      offerOnly = offerOnly, offertol = offertol, reltol = reltol, ...)

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
                     y, a = a, regr = regr, fn = logLikUlt, gr = logLikGradUlt ,
                     fixed = fvec, method = method, minOffer = minOffer, maxOffer = maxOffer,
                     offerOnly = offerOnly, offertol = offertol, reltol =
                     reltol, ...)
    }

    ## create a 'game' object to store output, with some extra attributes that
    ## are specific to the ultimatum model
    ans <- list()
    ans$coefficients <- results$estimate
    ans$vcov <- getGameVcov(results$hessian, fvec)
    ans$log.likelihood <-
        logLikUlt(results$estimate, y = y, acc = a, regr = regr, minOffer = minOffer,
        			maxOffer = maxOffer, offertol = offertol, offerOnly = offerOnly)
    ans$call <- cl
    ans$convergence <- list(method = method, iter = nIter(results), code =
                            results$code, message = results$message, gradient =
                            TRUE)
    ans$formulas <- formulas
    ans$link <- "logit"
    ans$type <- "private"
    ans$model <- mf
    ans$xlevels <- .getXlevels(attr(mf, "terms"), mf)
    ans$y <- ya
    ans$equations <- c("R1", "R2", "log(s1)", "log(s2)")
    attr(ans$equations, "hasColon") <- c(TRUE, TRUE, FALSE, FALSE)
    names(attr(ans$equations, "hasColon")) <- ans$equations
    ans$fixed <- fvec
    if (boot > 0)
        ans$boot.matrix <- bootMatrix
    ans$localID <- lid
    ans$outcome <- outcome
    ans$minOffer <- minOffer
    ans$maxOffer <- maxOffer
    ans$offertol <- offertol

    class(ans) <- c("game", "ultimatum")

    return(ans)
}
