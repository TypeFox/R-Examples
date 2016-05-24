`modifiedScoreStatistic` <-
function (fm, X, dispersion = 1) 
{
    y <- fm$y
    wt <- fm$prior
    LP <- fm$linear
    family <- fm$family
    probs <- family$linkinv(LP)
    dmu.deta <- family$mu.eta
    variance <- family$variance
    we <- c(wt * dmu.deta(LP)^2/variance(probs))
    W.X <- sqrt(we) * X
    XWXinv <- chol2inv(chol(crossprod(W.X)))
    hats <- diag(X %*% XWXinv %*% t(we * X))
    cur.model <- modifications(family, pl = fm$pl)(probs)
    mod.wt <- wt + c(hats * cur.model$at)
    y.adj <- (y * wt + hats * cur.model$ar)/mod.wt
    s.star <- t(c(dmu.deta(LP)/variance(probs)) * X) %*% ((y.adj - 
        probs) * mod.wt)
    t(s.star) %*% XWXinv %*% s.star
}
`penalizedDeviance` <-
function (fm, X, dispersion = 1) 
{
    Y <- fm$y
    LP <- fm$linear.predictor
    fam <- fm$family
    wt <- fm$prior.weights
    mu <- fm$fitted.values
    we <- fm$weights
    W.X <- sqrt(we) * X
    (sum(fam$dev.resid(Y, mu, wt)) - log(det(crossprod(W.X))))/dispersion
}
