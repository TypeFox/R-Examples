`RaoScoreStatistic` <-
function (fm, X, dispersion = 1) 
{
    fam <- fm$family
    mus <- fm$fitted
    tots <- fm$prior.weights
    variances <- tots * fam$variance(mus)
    dmu.deta <- tots * fam$mu.eta(fm$linear.predictor)
    Info <- crossprod((D <- X * dmu.deta)/sqrt(variances))/dispersion
    qScores <- t(D/variances) %*% ((fm$y - mus) * tots)/dispersion
    t(qScores) %*% chol2inv(chol(Info)) %*% qScores
}
`ordinaryDeviance` <-
function (fm, dispersion = 1) 
{
    LP <- fm$linear.predictor
    y <- fm$y
    fam <- fm$family
    mu <- fam$linkinv(LP)
    wt <- fm$prior.weights
    sum(fam$dev.resid(y, mu, wt))/dispersion
}
