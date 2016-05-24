## last modified August 29, 2002

mix <- function(mixdat, mixpar, dist = "norm", constr = list(conpi = "NONE", 
    conmu = "NONE", consigma = "NONE", fixpi = NULL, fixmu = NULL, 
    fixsigma = NULL, cov = NULL, size = NULL), emsteps = 1, usecondit = FALSE, 
    exptol = 5e-06, print.level = 0, ...) 
{
    k <- nrow(mixpar)
    if (usecondit & ncol(mixdat) - 2 != k) 
        stop("Conditional data are not consistent with mixpar.")
    if (usecondit & k == 1) 
        warning("Only 1 component, conditional data will not be used")
    if (is.na(match(dist, c("norm", "lnorm", "gamma", "weibull", 
        "binom", "nbinom", "pois")))) 
        stop(paste("Unknown distribution ", dist, ".", sep = ""))
    if (!testconstr(mixdat, mixpar, dist, constr)) 
        stop("Invalid constraints.")
    if (!testpar(mixpar, dist, constr)) 
        stop("Invalid parameters.")
    fitpar <- mixpar
    if (print.level > 0) 
        print(fitpar)
    estep <- function(emixdat, emixpar, edist, econstr) {
        prj.i <- grpintprob(emixdat, emixpar, edist, econstr)
        prji <- sweep(prj.i, 2, emixpar[, 1], "*")
        pri.j <- sweep(prji, 1, apply(prji, 1, sum), "/")
        pseudodat <- sweep(pri.j, 1, emixdat[, 2], "*")
        nplusi <- apply(pseudodat, 2, sum)
        pi <- nplusi/sum(nplusi)
        if (econstr$conpi == "PFX") {
            pi[econstr$fixpi] <- emixpar[econstr$fixpi, 1]
            pi[!econstr$fixpi] <- pi[!econstr$fixpi] * (1 - sum(pi[econstr$fixpi]))/sum(pi[!econstr$fixpi])
        }
        list(pi = pi, pseudomixdat = data.frame(x = emixdat[, 
            1], pseudodat = pseudodat))
    }
    mstep <- function(mspseudomixdat, msmixpar, msdist, msconstr, 
        msprint.level, ...) {
        multineg2llg <- function(mulpseudomixdat, mulmixpar, 
            muldist, mulconstr, multheta) {
            muly <- as.matrix(mulpseudomixdat[, -1])
            muln <- apply(muly, 2, sum)
            mulpara <- mixtheta2par(multheta, mulmixpar, mulconstr, 
                mixprop = FALSE)
            upmixpar <- data.frame(mulpi = mulmixpar[, 1], mulmu = mulpara[, 
                2], mulsigma = mulpara[, 3])
            mulparvalid <- testpar(upmixpar, muldist, mulconstr)
            if (mulparvalid) {
                mulp <- grpintprob(mulpseudomixdat, upmixpar, 
                  muldist, mulconstr)
                mulnp <- sweep(mulp, 2, muln, "*")
                mulyu <- as.vector(muly)
                mulnpu <- as.vector(mulnp)[mulyu > 0]
                mulyu <- mulyu[mulyu > 0]
                min(-2 * sum(mulyu * log(mulnpu/mulyu)), 1e+16)
            }
            else 1e+16
        }
        theta0 <- mixpar2theta(msmixpar, msconstr, mixprop = FALSE)
        msnlmobj <- nlm(multineg2llg, mulpseudomixdat = mspseudomixdat, 
            mulmixpar = msmixpar, muldist = msdist, mulconstr = msconstr, 
            p = theta0, print.level = msprint.level, ...)
        msk <- nrow(msmixpar)
        mspara <- mixtheta2par(msnlmobj$est, msmixpar, msconstr, 
            mixprop = FALSE)
        msmu <- mspara[, 2]
        mssigma <- mspara[, 3]
        data.frame(pi = msmixpar[, 1], mu = msmu, sigma = mssigma)
    }
    if (!usecondit & emsteps > 0 & k > 1 & (constr$conpi == "NONE" | 
        (constr$conpi == "PFX" & sum(constr$fixpi) < k - 1))) {
        for (i in 1:emsteps) {
            fitpi <- estep(mixdat, fitpar, dist, constr)
            fitpar[, 1] <- fitpi$pi
            fitpar <- mstep(fitpi$pseudomixdat, fitpar, dist, 
                constr, print.level, ...)
            if (print.level > 0) 
                print(fitpar)
        }
    }
    mixtheta0 <- mixpar2theta(fitpar, constr)
    mixlike <- function(lmixdat, lmixpar, ldist, lconstr, lmixtheta, 
        lusecondit, lexptol) {
        lk <- nrow(lmixpar)
        lm <- nrow(lmixdat)
        uppar <- mixtheta2par(lmixtheta, lmixpar, lconstr)
        parvalid <- testpar(uppar, ldist, lconstr)
        if (!parvalid) {
            like <- 1e+16
            ldf <- lm - 1
        }
        else {
            pmat <- grpintprob(lmixdat, uppar, ldist, lconstr)
            ln <- sum(lmixdat[, 2])
            joint <- sweep(pmat, 2, uppar[, 1], "*")
            mixed <- apply(joint, 1, sum)
            ly <- as.vector(as.matrix(lmixdat[, 2]))
            if (lusecondit & lk > 1) {
                conditional <- sweep(joint, 1, mixed, "/")
                conditdat <- as.matrix(lmixdat[, -(1:2)])
                subsam <- apply(conditdat, 1, sum)
                nconditpr <- sweep(conditional, 1, subsam, "*")
                conditdatu <- as.vector(conditdat)
                nconditpru <- as.vector(nconditpr)[conditdatu > 
                  0]
                conditdatu <- conditdatu[conditdatu > 0]
            }
            pu <- as.vector(mixed)[ly > 0]
            npu <- ln * pu
            yu <- ly[ly > 0]
            if (lusecondit) {
                like <- min(-2 * (sum(yu * log(npu/yu)) + sum(conditdatu * 
                  log(nconditpru/conditdatu))), 1e+16)
                ldf <- sum(mixed > lexptol) - 1 + sum(nconditpr > 
                  lexptol) - sum(subsam > 0) - length(lmixtheta)
            }
            else {
                like <- min(-2 * sum(yu * log(npu/yu)), 1e+16)
                ldf <- sum(mixed > lexptol) - 1 - length(lmixtheta)
            }
        }
        attr(like, "df") <- ldf
        like
    }
    nlmobj <- nlm(mixlike, lmixdat = mixdat, lmixpar = fitpar, 
        ldist = dist, lconstr = constr, lusecondit = usecondit, 
        lexptol = exptol, p = mixtheta0, hessian = TRUE, print.level = print.level, 
        ...)
    if (nlmobj$code == 3) 
        warning("The optimization process terminated because either the estimates are approximate local optimal solution or steptol is too small")
    else if (nlmobj$code == 4) 
        warning("The optimization process terminated because iteration limit exceeded")
    else if (nlmobj$code == 5) 
        warning("The optimization process terminated because either the function is unbounded below or stepmax is too small")
    if (nlmobj$minimum == 1e+16) 
        warning("Invalid parameters")
    param <- mixtheta2par(nlmobj$est, mixpar, constr)
    covmat <- function(covmixpar, covconstr, covhessian) {
        invmat <- try(solve(covhessian/2))
        if (inherits(invmat, "try-error")) 
            invmat <- matrix(NA, nrow = nrow(covhessian), ncol = ncol(covhessian))
        dr <- nrow(invmat)
        covk <- nrow(covmixpar)
        if (covconstr$conpi == "NONE") 
            lpi <- covk - 1
        else if (covconstr$conpi == "PFX" & sum(covconstr$fixpi) < 
            covk - 1) {
            pi.e <- covmixpar[!covconstr$fixpi, 1]
            lpi <- length(pi.e) - 1
        }
        else lpi <- 0
        if (covconstr$conmu == "NONE") 
            mu.e <- covmixpar[, 2]
        else if (covconstr$conmu == "MFX") 
            mu.e <- covmixpar[!covconstr$fixmu, 2]
        else if (covconstr$conmu == "MEQ") 
            mu.e <- covmixpar[1, 2]
        else if (covconstr$conmu == "MES") 
            mu.e <- covmixpar[1:2, 2]
        else if (covconstr$conmu == "MGC") 
            mu.e <- covmixpar[1:3, 2]
        lmu <- length(mu.e)
        if (covconstr$consigma == "NONE") 
            sigma.e <- covmixpar[, 3]
        else if (covconstr$consigma == "SFX") 
            sigma.e <- covmixpar[!constr$fixsigma, 3]
        else if (!is.na(match(covconstr$consigma, c("SEQ", "CCV")))) 
            sigma.e <- covmixpar[1, 3]
        else if (!is.na(match(covconstr$consigma, c("FCV", "BINOM", 
            "NBINOM", "POIS")))) 
            sigma.e <- NULL
        lsigma <- length(sigma.e)
        sigmat <- diag(c(rep(1, dr - lsigma), sigma.e))
        vmat <- sigmat %*% invmat %*% sigmat
        se <- sqrt(diag(vmat))
        pi.se <- rep(NA, covk)
        if (lpi > 0) {
            if (covconstr$conpi == "PFX") 
                pi.se[!covconstr$fixpi] <- c(se[1:lpi], sqrt(sum(vmat[1:lpi, 
                  1:lpi])))
            else pi.se[1:covk] <- c(se[1:lpi], sqrt(sum(vmat[1:lpi, 
                1:lpi])))
        }
        mu.se <- rep(NA, covk)
        if (lmu > 0) {
            if (covconstr$conmu == "MFX") 
                mu.se[!covconstr$fixmu] <- se[(lpi + 1):(lpi + 
                  lmu)]
            else mu.se[1:lmu] <- se[(lpi + 1):(lpi + lmu)]
        }
        sigma.se <- rep(NA, covk)
        if (lsigma > 0) {
            if (covconstr$consigma == "SFX") 
                sigma.se[!covconstr$fixsigma] <- se[(lpi + lmu + 
                  1):(lpi + lmu + lsigma)]
            else sigma.se[1:lsigma] <- se[(lpi + lmu + 1):(lpi + 
                lmu + lsigma)]
        }
        list(vmat = vmat, mixse = data.frame(pi.se = pi.se, mu.se = mu.se, 
            sigma.se = sigma.se))
    }
    vmat <- covmat(param, constr, nlmobj$hessian)
    chisq <- nlmobj$minimum
    df <- attr(mixlike(mixdat, fitpar, dist, constr, nlmobj$est, 
        usecondit, exptol), "df")
    P <- 1 - pchisq(chisq, df)
    mixfit <- list(parameters = param, se = vmat$mixse, distribution = dist, 
        constraint = constr, chisq = chisq, df = df, P = P, vmat = vmat$vmat, 
        mixdata = mixdat, usecondit = usecondit)
    class(mixfit) <- "mix"
    mixfit
}
