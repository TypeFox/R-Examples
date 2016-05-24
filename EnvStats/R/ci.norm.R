ci.norm <-
function (muhat, sdhat, n, ci.type, alpha) 
{
    sd.muhat <- sdhat/sqrt(n)
    df <- n - 1
    ret.obj <- ci.normal.approx(muhat, sd.muhat, n, df, ci.type, 
        alpha)
    ret.obj$parameter <- "mean"
    ret.obj$method <- "Exact"
    ret.obj
}
