zlm <-
function (formula, data = NULL, subset = NULL, g = "UIP") 
{
    thiscall = match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if (is.matrix(formula)) {
        mf <- model.frame(as.data.frame(formula, drop.unused.levels = TRUE))
    }
    else {
        mf <- eval(mf, parent.frame())
    }
    yXdata = as.matrix(mf)
    N = nrow(yXdata)
    K = ncol(yXdata) - 1
    dmdata = yXdata - matrix(colMeans(yXdata), N, K + 1, byrow = TRUE)
    yty = c(crossprod(dmdata[, 1]))
    olsres = .ols.terms2(positions = rep(TRUE, K), yty = yty, 
        N = N, K = K, XtX.big = crossprod(dmdata[, -1, drop = FALSE]), 
        Xty.big = c(crossprod(dmdata[, -1, drop = FALSE], dmdata[, 
            1])))$full.results()
    if (is.list(g)) {
        if (any(is.element(names(g), "gtype"))) 
            gprior.info = g
        else stop("Please provide a proper g-prior. see help(zlm)")
    }
    gprior.info = .choose.gprior(g = g, N = N, K = K, return.g.stats = TRUE, 
        yty = yty)
    lprobcalc = gprior.info$lprobcalc
    zres = lprobcalc$lprob.all(ymy = olsres$ymy, k = K, bhat = olsres$bhat, 
        diag.inverse = olsres$diag.inverse)
    betas = c(zres$b1)
    betas2 = c(zres$b2)
    alpha = mean(yXdata[, 1]) - c(crossprod(betas, colMeans(yXdata)[-1]))
    fitval = c(yXdata[, -1, drop = FALSE] %*% betas) + alpha
    resids = yXdata[, 1] - fitval
    if (gprior.info$is.constant) {
        gprior.info$shrinkage.moments = 1 - 1/(1 + gprior.info$g)
    }
    else {
        gprior.info$shrinkage.moments = zres$otherstats
    }
    mt <- attr(mf, "terms")
    alphabeta = c(alpha, betas)
    names(alphabeta) <- c("(Intercept)", attr(mt, "term.labels"))
    res = list()
    res$coefficients <- alphabeta
    res$residuals <- resids
    res$rank <- K + 1
    res$fitted.values <- fitval
    res$df.residual <- N - K - 1
    res$xlevels <- .getXlevels(mt, mf)
    res$call <- thiscall
    res$terms <- mt
    res$model <- mf
    res$na.action <- attr(mf, "na.action")
    res$coef2moments <- c(NA, betas2)
    res$marg.lik <- zres$lprob
    res$gprior.info <- gprior.info
    res$olsres <- olsres
    res$zres <- zres
    class(res) = c("zlm", "lm")
    return(res)
}
