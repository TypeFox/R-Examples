###############################################################################
## ordinalNet package includes:
###############################################################################
## getObjfun
## getIC
## invLogit
## softThresh
## yFactorToMatrix
## makeOrdinalXLS
## makeOrdinalLink
## makeConvFun
## mirlsNetSV
## mirlsNet
## ordinalNet
## print.mirlsNetFit
## coef.mirlsNetFit
## predict.mirlsNetFit
## print.ordinalNetFit
## coef.ordinalNetFit
## predict.ordinalNetFit
###############################################################################

getObjfun <- function(xLS, yMat, linkfun, betaHat, alpha, lambda, penalizeID)
{
    nObs <- sum(yMat)
    etaLS <- lapply(xLS, function(x) x %*% betaHat)
    pLS <- lapply(etaLS, function(eta) linkfun$h(eta))
    pLS <- lapply(pLS, function(p) c(p, 1-sum(p)))
    pMat <- do.call(rbind, pLS)
    if (any(pMat<0))
    {
        loglik <- -Inf
    } else
    {
        llMat <- yMat * log(pMat)
        llMat[yMat==0] <- 0  # -Inf*0 = 0
        loglik <- sum(llMat)
    }

    l1 <- sum(abs(betaHat)*penalizeID)
    l22 <- sum(betaHat^2*penalizeID)
    pen1 <- ifelse(l1==0, 0, alpha*lambda*l1)
    pen2 <- ifelse(l22==0, 0, .5*(1-alpha)*l22)
    objfun <- -1/nObs*loglik + pen1 + pen2
    objfun
}

getIC <- function(xLS, yMat, linkfun, betaHat)
{
    nObs <- sum(yMat)
    etaLS <- lapply(xLS, function(x) x %*% betaHat)
    pLS <- lapply(etaLS, function(eta) linkfun$h(eta))
    pLS <- lapply(pLS, function(p) c(p, 1-sum(p)))
    pMat <- do.call(rbind, pLS)
    if (any(pMat<0))
    {
        loglik <- -Inf
    } else
    {
        llMat <- yMat * log(pMat)
        llMat[yMat==0] <- 0  # -Inf*0 = 0
        loglik <- sum(llMat)
    }

    df <- sum(betaHat != 0)
    aic <- 2*df - 2*loglik
    bic <- log(nObs)*df - 2*loglik

    y <- colSums(yMat)
    pHatNull <-  y / nObs
    loglikNull <- crossprod(y, log(pHatNull))
    devPct <- 1 - loglik / loglikNull

    c(df=df, aic=aic, bic=bic, devPct=devPct)
}

invLogit <- function(x) 1/(1+exp(-x))

softThresh <- function(z, g) sign(z)*max(0, abs(z)-g)

yFactorToMatrix <- function(y)
{
    nObs <- length(y)
    nLev <- length(levels(y))
    yMat <- matrix(0, nrow=nObs, ncol=nLev, dimnames=list(NULL, levels(y)))
    yInt <- as.integer(y)
    yMat[cbind(1:nObs, yInt)] <- 1
    yMat
}

makeOrdinalXLS <- function(x, nLev)
{
    n <- nrow(x)
    d <- diag(nLev-1)
    xLS <- lapply(1:n, function(i) cbind(d, x[rep(i, nLev-1),,drop=FALSE]))
    xLS
}

#' @importFrom stats make.link
makeOrdinalLink <- function(link)
{
    oklinks <- c("logit", "probit", "cloglog", "cauchit")
    if (!link %in% oklinks)
        stop("Invalid link specified for ordinal multinomial regression model")
    linkfun <- make.link(link)

    ## g(p) is the multivariate link function
    g <- function(p)
    {
        gTemp <- linkfun$linkfun
        gamma <- cumsum(p)
        g <- gTemp(gamma)
        g
    }

    ## h(eta) is the multivariate inverse link function
    h <- function(eta)
    {
        hTemp <- linkfun$linkinv
        hEta <- hTemp(eta)
        k <- length(eta)
        h <- c(hEta[1], hEta[-1] - hEta[-k])
        h
    }

    ## getQ(p) returns the Jacobian of the inverse link function at p
    getQ <- function(eta)
    {
        k <- length(eta)
        muEta <- linkfun$mu.eta
        ## nrow=k makes diagQ a matrix in the case when k=1
        diagQ <- diag(muEta(eta), nrow=k)
        offDiagQ <- rbind(rep(0, k), -diagQ[-k,])
        q <- diagQ + offDiagQ
        q
    }

    list(g=g, h=h, getQ=getQ)
}

## norm used to check convergence of beta and gamma
makeConvFun <- function(convNorm)
{
    convFun <- if (convNorm==Inf) {
        function(x) max(abs(x))
    } else {
        function(x) mean(abs(x)^convNorm)^(1/convNorm)
    }
    convFun
}

## Main function - can be used for non-ordinal models
## yMat used for objective function; yLS used for algorithm
mirlsNetSV <- function(xLS, yMat, yLS, lambda, alpha=1, penalizeID, positiveID=NULL, linkfun, betaStart,
                       quadApp=NULL, trace=FALSE, epsOut=1e-3, epsIn=1e-3, maxiterOut=Inf,
                       maxiterIn=Inf, pMin=1e-8, betaMin=1e-8, convNorm=Inf, ordinalCheck=FALSE)
{
    if (is.null(positiveID)) positiveID <- rep(FALSE, ncol(xLS[[1]]))

    nObs <- sum(unlist(yLS))
    nLev <- length(yLS[[1]])
    nVar <- ncol(xLS[[1]])
    convFun <- makeConvFun(convNorm)

    betaHat <- betaStart
    objfun <- getObjfun(xLS, yMat, linkfun, betaHat, alpha, lambda, penalizeID)
    convOut <- FALSE
    iterOut <- 0
    while (!convOut && iterOut < maxiterOut)
    {
        iterOut <- iterOut + 1
        betaHatOldOut <- betaHat

        if (iterOut==1 && !is.null(quadApp))
        {
            activeID <- quadApp$activeID
            xwx <- quadApp$xwx
            term <- quadApp$term
            betaHat <- quadApp$betaHat
            xwxActive <- quadApp$xwxActive
            termActive <- quadApp$termActive
            betaHatActive <- quadApp$betaHatActive
            penalizeIDActive <- quadApp$penalizeIDActive
            positiveIDActive <- quadApp$positiveIDActive
        } else
        {
            ## Update Quadratic Approximation:
            #######################################################################
            ## Note: yLS, etaLS, cdfLS, pLS are lists of vectors
            ## yLS vectors have length nLev while the other vectors have length nLev-1
            etaLS <- lapply(xLS, function(x) as.vector(x %*% betaHat))
            pLS <- lapply(etaLS, function(eta) linkfun$h(eta))
            pkPlusOneLS <- lapply(pLS, function(p) 1-sum(p))
            pHatMin <- min(unlist(pLS), unlist(pkPlusOneLS))
            if (pHatMin < pMin) break
            cdfLS <- lapply(pLS, function(p) cumsum(p))
            dLS <- mapply(y=yLS, p=pkPlusOneLS, function(y, p) y[nLev] / p, SIMPLIFY=FALSE)
            uLS <- mapply(y=yLS, p=pLS, function(y, p) y[-nLev] / p, SIMPLIFY=FALSE)
            qLS <- lapply(etaLS, linkfun$getQ)
            jj <- matrix(1, ncol=nLev-1, nrow=nLev-1)
            sigInvLS <- mapply(y=yLS, p=pLS, function(y, p) sum(y) * (diag(1/p, nrow=nLev-1) + jj / (1-sum(p))), SIMPLIFY=FALSE)
            wLS <- mapply(q=qLS, s=sigInvLS, function(q, s) crossprod(q, s) %*% q, SIMPLIFY=FALSE)

            xwxLS <- mapply(x=xLS, w=wLS, function(x, w) crossprod(x, w) %*% x, SIMPLIFY=FALSE)
            xqudLS <- mapply(x=xLS, q=qLS, u=uLS, d=dLS, function(x, q, u, d) crossprod(x, crossprod(q, u-d)), SIMPLIFY=FALSE)
            xwx <- Reduce("+", xwxLS)
            xqud <- Reduce("+", xqudLS)

            ## Update Active Set
            #######################################################################
            activeID <- diag(xwx)!=0 & (!penalizeID | betaHat!=0)
            betaHat[!activeID] <- 0
            betaHatActive <- betaHat[activeID]
            term <- xwx[, activeID, drop=FALSE] %*% betaHatActive + xqud
            xwxActive <- xwx[activeID, activeID, drop=FALSE]
            termActive <- term[activeID]
            penalizeIDActive <- penalizeID[activeID]
            positiveIDActive <- positiveID[activeID]
        }

        ## Coordinate Descent inner loop:
        #######################################################################
        convKKT <- FALSE
        while (!convKKT)
        {
            iterIn <- 0
            convIn <- FALSE
            while (!convIn && iterIn<maxiterIn)
            {
                iterIn <- iterIn + 1
                betaHatOldIn <- betaHat
                for (i in 1:sum(activeID))
                {
                    num <- termActive[i] - xwxActive[i, -i, drop=FALSE] %*% betaHatActive[-i]
                    den <- xwxActive[i, i]
                    if (penalizeIDActive[i])
                    {
                        betaHatActive[i] <- softThresh(num, nObs*lambda*alpha) / (den+nObs*lambda*(1-alpha))
                    } else
                    {
                        betaHatActive[i] <- num / den
                    }
                    if (positiveIDActive[i]) betaHatActive[i] <- max(0, betaHatActive[i])
                }
                betaHat[activeID] <- betaHatActive
                betaHat[abs(betaHat) < betaMin] <- 0  # for terms converging slowly
                relChangeIn <- (betaHat - betaHatOldIn) / betaHatOldIn
                relChangeIn[betaHat == betaHatOldIn] <- 0  # prevents div/0
                difIn <- convFun(relChangeIn)
                convIn <- difIn < epsIn
            }
            kkt <- rep(TRUE, nVar)
            kktInactive <- term[!activeID] - xwx[!activeID, activeID, drop=FALSE] %*% betaHatActive
            kktInactive[!positiveID[!activeID]] <- abs(kktInactive[!positiveID[!activeID]])
            kkt[!activeID] <- kktInactive <= nObs*lambda*alpha
            convKKT <- all(kkt)
            if (!convKKT)
            {
                activeID[!kkt] <- TRUE
                xwxActive <- xwx[activeID, activeID, drop=FALSE]
                termActive <- term[activeID]
                betaHatActive <- betaHat[activeID]
                penalizeIDActive <- penalizeID[activeID]
                positiveIDActive <- positiveID[activeID]
            }
        }
        #######################################################################
        if (!convIn) break

        ## Checks
        ########################################################################
        # Intercepts should be non-decreasing (Only applies to ordinal models)
        if (ordinalCheck) betaHat[1:(nLev-1)] <- sort(betaHat[1:(nLev-1)])

        # If loglik does not improve, take half step size
        # Note: Do not loop until ll>llOld. This can get stuck.
        objfunOld <- objfun
        objfun <- getObjfun(xLS, yMat, linkfun, betaHat, alpha, lambda, penalizeID)
        if (objfun > objfunOld) betaHat <- (betaHat + betaHatOldOut) / 2
        ########################################################################

        relChangeOut <- (betaHat - betaHatOldOut) / betaHatOldOut
        relChangeOut[betaHat == betaHatOldOut] <- 0  # prevents div/0
        difOut <- convFun(relChangeOut)
        convOut <- difOut < epsOut

        if (trace) cat("Outer Iteration ", iterOut, ":  ", iterIn, " inner iterations,  ",
                       signif(difOut, 2), " relative change\n", sep='')
    }
    quadApp <- list(activeID=activeID, xwx=xwx, term=term, betaHat=betaHat,
                    xwxActive=xwxActive, termActive=termActive, betaHatActive=betaHatActive,
                    penalizeIDActive=penalizeIDActive, positiveIDActive=positiveIDActive)
    mirlsNetSV <- list(converged=convOut, iterOut=iterOut, betaHat=betaHat, quadApp=quadApp)
    mirlsNetSV
}

## Wrapper for mirlsNetSV
mirlsNet <-function(xLS, yMat, alpha=1, penalizeID, positiveID=rep(FALSE, nVar),
                    linkfun, betaStart, lambdaVals=NULL, nLambda=100,
                    lambdaMinRatio=ifelse(nObs<nVar, 1e-2, 1e-4), alphaMin=0.01,
                    trace=FALSE, epsOut=1e-3, epsIn=1e-3, maxiterOut=Inf,
                    maxiterIn=Inf, pMin=1e-8, betaMin=1e-8, convNorm=Inf, ordinalCheck=FALSE)
{
    nObs <- sum(yMat)
    nVar <- ncol(xLS[[1]])
    if (!is.null(lambdaVals)) lambdaVals <- sort(lambdaVals, decreasing=TRUE)
    arglist <- as.list(environment())
    arglist[c("nObs", "nVar")] <- NULL

    yLS <- lapply(1:nrow(yMat), function(i) yMat[i,])
    mods <- vector("list", length=ifelse(!is.null(lambdaVals), length(lambdaVals), nLambda))
    quadApp <- NULL
    if (is.null(lambdaVals))
    {
        if (trace) cat("Lambda", 1, " of ", nLambda, "\n", sep='')
        mods[[1]] <- mirlsNetSV(xLS, yMat, yLS, lambda=Inf, alpha, penalizeID, positiveID,
                                linkfun, betaStart, quadApp, trace, epsOut, epsIn, maxiterOut,
                                maxiterIn, pMin, betaMin, convNorm, ordinalCheck)
        betaStart <- mods[[1]]$betaHat
        quadApp <- mods[[1]]$quadApp
        kktInactive <- with(quadApp, term[penalizeID] - xwx[penalizeID, activeID, drop=FALSE] %*% betaHatActive)
        kktInactive[!positiveID[penalizeID]] <- abs(kktInactive[!positiveID[penalizeID]])
        kktInactive <- kktInactive / (nObs*max(alpha, alphaMin))
        lambdaMax <- max(kktInactive, 0)
        lambdaMin <- lambdaMax * lambdaMinRatio
        lambdaVals <- exp(seq(log(lambdaMax), log(lambdaMin), length.out=nLambda))
    }
    for (i in which(sapply(mods, is.null)))
    {
        if (trace) cat("Lambda", i, " of ", length(lambdaVals), "\n")
        lambda <- lambdaVals[i]
        mods[[i]] <- mirlsNetSV(xLS, yMat, yLS, lambda, alpha, penalizeID, positiveID,
                                linkfun, betaStart, quadApp, trace, epsOut, epsIn, maxiterOut,
                                maxiterIn, pMin, betaMin, convNorm, ordinalCheck)
        betaStart <- mods[[i]]$betaHat
        quadApp <- mods[[i]]$quadApp
        if (!mods[[i]]$converged)
        {
            if (i==1) stop("Coefficient estimates diverging. Try larger lambdaVals or set lambdaVals=NULL to automatically generate lambda sequence.")
            break
        }
    }
    mods <- mods[1:(i-!mods[[i]]$converged)]
    lambdaVals <- lambdaVals[seq_along(mods)]
    arglist$lambdaVals <- lambdaVals
    betaHat <- do.call(rbind, lapply(mods, function(m) m$betaHat))
    icVals <- do.call(rbind, lapply(mods, function(m) getIC(xLS, yMat, linkfun, m$betaHat)))
    converged <- do.call(c, lapply(mods, function(m) m$converged))
    iterOut <- do.call(c, lapply(mods, function(m) m$iterOut))
    mirlsNetFit <- c(arglist, list(iterOut=iterOut, converged=converged, icVals=icVals, betaHat=betaHat))
    class(mirlsNetFit) <- "mirlsNetFit"
    mirlsNetFit
}

## ordinalNet is a wrapper for mirlsNet
#' Penalized ordinal regression
#'
#' Fits ordinal regression models with elastic net penalty.
#' Supported models include cumulative logit, probit, cauchit, and complementary
#' log-log. The regularization path is computed at a grid of values for the
#' regularizaton parameter \code{lambda}. The algorithm uses Fisher Scoring
#' with Coordinate Descent updates.
#'
#' @param x Covariate matrix
#' @param y Response variable, either an ordered factor or a matrix
#' where each row is a multinomial vector of counts
#' @param alpha The elastic net mixing parameter, with \eqn{0\le\alpha\le1}.
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} is the ridge penalty.
#' See "Details".
#' @param standardize If \code{standardize=TRUE}, the predictor variables are
#' scaled to have unit variance. Coefficient estimates are returned on the
#' original scale.
#' @param penalizeID Logical vector indicating whether each coefficient should
#' be penalized. Default is TRUE for all coefficients.
#' @param positiveID Logical vector indicating whether each coefficient should
#' be constrained to be non-negative. Default is FALSE for all coefficients.
#' @param link Specifies the link function. The options supported are logit,
#' probit, complementary log-log, and cauchit.
#' @param lambdaVals An optional user-specified lambda sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nLambda}
#' and \code{lambdaMinRatio}.
#' @param nLambda The number of \code{lambda} values in the solution path.
#' Default is 100.
#' @param lambdaMinRatio Smallest value for \code{lambda} as a fraction of the maximum
#' lambda. (The maximum lambda is the smallest value that sets all penalized
#' coefficients to zero.)
#' @param alphaMin If \code{alpha < alphaMin}, then \code{alphaMin} is used
#' to calculate the maximum lambda value.
#' @param trace If \code{trace=TRUE} the algorithm progress is printed to the terminal.
#' @param epsOut Convergence threshold for the algorithm's outer loop. The outer
#' loop optimizes optimizes the Fisher Scoring quadratic approximation to the
#' penalized log likelihood.
#' @param epsIn Convergence threshold for the algorithm's inner loop. The inner
#' loop cycles through and updates the coefficient estimates using coordinate
#' descent. Each cycle is one iteration.
#' @param maxiterOut Maximum number of outer loop iterations.
#' @param maxiterIn Maximum number of inner loop iterations.
#' @param pMin If for any observation, the fitted probability for a response
#' category falls below pMin, the algorithm is terminated. This can occur for
#' small \code{lambda} values as the coefficient estimates diverge to \eqn{+/-\infty}.
#' @param betaMin If a coefficient estimate falls below \code{betaMin}, it is
#' set to zero. This improves the stability and speed of the fitting algorithm.
#' @param convNorm The Lp norm of the coefficient estimate relative changes is
#' computed after each iteration of the inner or outer loop. Convergence of the
#' loop is assessed by comparing this value to the convergence threshold
#' (\code{epsIn} or \code{epsOut}). \code{convNorm=Inf} is the default and
#' represents the L-\eqn{\infty} norm (maximum relative change).
#' @param zetaStart Optional user-specified starting values for the intercept
#' terms (must be non-decreasing). Default is a uniform sequency from -1 to 1.
#' @param thetaStart Optional user-specified starting values for the non-intercept
#' terms. Default is zero for all coefficients.
#' @return An object with S3 class "ordinalNetFit".
#' @details
#' The ordinal model has the form
#' \deqn{g(P(y\le j|x)) = Intercept[j] + x\beta,}
#' where \eqn{g(.)} is a link function, most commonly logit.
#'
#' The elastic net penalty is defined as
#' \deqn{\lambda{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1}.}
#'
#' The objective function is
#' \deqn{-1/N*loglik + penalty.}
#' @examples
#' set.seed(10)
#' x <- matrix(rnorm(50*5), ncol=5)
#' beta <- c(1, 0, 0, 0, 0)
#' intercepts <- c(-1, 1)
#' xb <- x %*% beta
#' eta <- cbind(xb+intercepts[1], xb+intercepts[2])
#' probMatCumul <- 1 / (1 + exp(-eta))
#' probMat <- cbind(probMatCumul, 1) - cbind(0, probMatCumul)
#' y <- apply(probMat, 1, function(p) sample(1:3, size=1, prob=p))
#' y <- as.factor(y)
#' fit <- ordinalNet(x, y)
#' print(fit)
#' coef(fit)
#' predict(fit, type="class")
#' predict(fit, type="prob")
#' @export
#' @importFrom stats sd
ordinalNet <- function(x, y, alpha=1, standardize=TRUE, penalizeID=NULL, positiveID=NULL,
                       link=c("logit", "probit", "cloglog", "cauchit"),
                       lambdaVals=NULL, nLambda=100, lambdaMinRatio=ifelse(nObs<nVar, 1e-2, 1e-4),
                       alphaMin=0.01, trace=FALSE, epsOut=1e-3, epsIn=1e-3,
                       maxiterOut=Inf, maxiterIn=Inf, pMin=1e-20, betaMin=1e-8,
                       convNorm=Inf, zetaStart=NULL, thetaStart=NULL)
{
    yMat <- if (is.factor(y))  yFactorToMatrix(y) else y
    nObs <- nrow(yMat)
    nVar <- ncol(x)
    nLev <- ncol(yMat)
    if (is.null(zetaStart)) zetaStart <- seq(-1, 1, length.out=nLev-1)
    if (is.null(thetaStart)) thetaStart <- rep(0, nVar)
    betaStart <- c(zetaStart, thetaStart)
    if (is.null(penalizeID)) penalizeID <- rep(TRUE, nVar)
    if (is.null(positiveID)) positiveID <- rep(FALSE, nVar)
    link <- match.arg(link)
    linkfun <- makeOrdinalLink(link)
    penalizeID <- c(rep(FALSE, nLev-1), penalizeID)
    positiveID <- c(rep(FALSE, nLev-1), positiveID)

    if (standardize)
    {
        xWts <- rowSums(yMat)
        xSD <- apply(x[rep(1:nrow(x), xWts), , drop=FALSE], 2, sd) * sqrt((nObs-1)/nObs)
        xSD[xSD==0] <- 1
        xScaled <- x / rep(xSD, each=nrow(x))
        xLS <- makeOrdinalXLS(xScaled, nLev)
    } else
    {
        xLS <- makeOrdinalXLS(x, nLev)
    }


    ordinalCheck <- TRUE
    mirlsNetFit <- mirlsNet(xLS, yMat, alpha, penalizeID, positiveID,
                            linkfun, betaStart, lambdaVals, nLambda,
                            lambdaMinRatio, alphaMin,
                            trace, epsOut, epsIn, maxiterOut,
                            maxiterIn, pMin, betaMin, convNorm, ordinalCheck)

    nonint <- mirlsNetFit$betaHat[ , -(1:nLev-1), drop=FALSE]
    int <- mirlsNetFit$betaHat[ , 1:(nLev-1), drop=FALSE]
    if (standardize) nonint <- nonint / rep(xSD, each=nrow(nonint))
    coef <- cbind(int, nonint)
    intNames <- paste0("Intercept", 1:(nLev-1))
    xNames <- if (is.null(colnames(x))) paste0("X", 1:nVar) else colnames(x)
    colnames(coef) <- c(intNames, xNames)

    ordinalNetFit <- list(coef=coef, mirlsNetFit=mirlsNetFit)
    class(ordinalNetFit) <- "ordinalNetFit"
    ordinalNetFit
}

print.mirlsNetFit <- function(x, ...) print(do.call(cbind, x[c("icVals", "lambdaVals")]))

coef.mirlsNetFit <- function(object, whichLambda=NULL, criteria=c("aic", "bic"), ...)
{
    criteria <- match.arg(criteria)
    if (is.null(whichLambda)) whichLambda <- which.min(object$icVals[,criteria])
    betaHat <- object$betaHat[whichLambda,]
    betaHat
}

predict.mirlsNetFit <- function(object, newx=NULL, whichLambda=NULL, criteria=c("aic", "bic"), type=c("class", "prob"), ...)
{
    type <- match.arg(type)
    criteria <- match.arg(criteria)
    if (is.null(whichLambda)) whichLambda <- which.min(object$icVals[,criteria])
    betaHat <- object$betaHat[whichLambda,]
    xLS <- if (is.null(newx)) object$xLS else makeOrdinalXLS(newx, ncol(object$yMat))
    etaLS <- lapply(xLS, function(x) x %*% betaHat)
    probLS <- lapply(etaLS, object$linkfun$h)
    probLS <- lapply(probLS, function(p) c(p, 1-sum(p)))
    probMat <- do.call(rbind, probLS)
    attr(probMat, "dimnames") <- NULL
    class <- apply(probMat, 1, which.max)
    switch(type, class=class, prob=probMat)
}

#' Print method for an "ordinalNetFit" object.
#'
#' Provides a model fit summary in matrix form. For each \code{lambda} value in
#' the solution path, the following information is included: degrees of freedom
#' (number of nonzero parameters), AIC, BIC, and percent deviance explained.
#'
#' @param x An "ordinalNetFit" S3 object.
#' @param ... Not used. Additional print arguments.
#' @return Silently returns the model fit summary.
#' @examples
#' set.seed(10)
#' x <- matrix(rnorm(50*5), ncol=5)
#' beta <- c(1, 0, 0, 0, 0)
#' intercepts <- c(-1, 1)
#' xb <- x %*% beta
#' eta <- cbind(xb+intercepts[1], xb+intercepts[2])
#' probMatCumul <- 1 / (1 + exp(-eta))
#' probMat <- cbind(probMatCumul, 1) - cbind(0, probMatCumul)
#' y <- apply(probMat, 1, function(p) sample(1:3, size=1, prob=p))
#' y <- as.factor(y)
#' fit <- ordinalNet(x, y)
#' print(fit)
#' @export
print.ordinalNetFit <- function(x, ...) print.mirlsNetFit(x$mirlsNetFit)

#' Extracts fitted coefficients from an "ordinalNetFit" object.
#'
#' @param object An "ordinalNetFit" S3 object.
#' @param whichLambda Optional index number(s) of the desired \code{lambda} within
#' the sequence of \code{lambdaVals}.
#' \code{whichLambda=1:nLambda} will return the entire solution path.
#' @param criteria Selects the best \code{lambda} value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param ... Not used. Additional coef arguments.
#' @return Vector of fitted coefficients.
#' @examples
#' set.seed(10)
#' x <- matrix(rnorm(50*5), ncol=5)
#' beta <- c(1, 0, 0, 0, 0)
#' intercepts <- c(-1, 1)
#' xb <- x %*% beta
#' eta <- cbind(xb+intercepts[1], xb+intercepts[2])
#' probMatCumul <- 1 / (1 + exp(-eta))
#' probMat <- cbind(probMatCumul, 1) - cbind(0, probMatCumul)
#' y <- apply(probMat, 1, function(p) sample(1:3, size=1, prob=p))
#' y <- as.factor(y)
#' fit <- ordinalNet(x, y)
#' coef(fit)
#' @export
coef.ordinalNetFit <- function(object, whichLambda=NULL, criteria=c("aic", "bic"), ...)
{
    object$mirlsNetFit$betaHat <- object$coef
    coef.mirlsNetFit(object$mirlsNetFit, whichLambda=whichLambda, criteria=criteria)
}

#' Predict method for an "ordinalNetFit" object
#'
#' Obtains predicted class numbers or fitted class probabilities for either
#' the model matrix used to fit the model or a new model matrix.
#'
#' @param object An "ordinalNetFit" S3 object.
#' @param newx Optional covariate matrix. If the model was fit with a centered
#' and scaled matrix, then \code{newx} should be centered and scaled as well.
#' @param whichLambda Optional index number of the desired \code{lambda} within
#' the sequence of \code{lambda} values in the solution path.
#' @param criteria Selects the best \code{lambda} value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param type Specifies whether to return a vector of predicted class numbers
#' or a matrix of fitted class probabilities.
#' @param ... Not used. Additional predict arguments.
#' @return A vector of predicted class numbers or a matrix of fitted class probabilities,
#' depending on \code{type}.
#' @examples
#' set.seed(10)
#' x <- matrix(rnorm(50*5), ncol=5)
#' beta <- c(1, 0, 0, 0, 0)
#' intercepts <- c(-1, 1)
#' xb <- x %*% beta
#' eta <- cbind(xb+intercepts[1], xb+intercepts[2])
#' probMatCumul <- 1 / (1 + exp(-eta))
#' probMat <- cbind(probMatCumul, 1) - cbind(0, probMatCumul)
#' y <- apply(probMat, 1, function(p) sample(1:3, size=1, prob=p))
#' y <- as.factor(y)
#' fit <- ordinalNet(x, y)
#' predict(fit, type="class")
#' predict(fit, type="prob")
#' @export
predict.ordinalNetFit <- function(object, newx=NULL, whichLambda=NULL, criteria=c("aic", "bic"), type=c("class", "prob"), ...)
{
    if (!is.null(newx)) object$mirlsNetFit$betaHat <- object$coef
    predict.mirlsNetFit(object$mirlsNetFit, newx=newx, whichLambda=whichLambda, criteria=criteria, type=type)
}
