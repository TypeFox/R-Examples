`orglellipse` <-
    function(object, groups, display = "sites", w = weights(object, display),
             kind = c("sd", "se"), conf, choices = 1:3, alpha = 0.3,
             col = "red", ...)
{
    weights.default <- function(object, ...) NULL
    kind <- match.arg(kind)
    x <- scores(object, display = display, choices = choices, ...)
    groups <- as.factor(groups)
    ## evaluate weights
    w <- eval(w)
    if (is.null(w) || length(w) == 1)
        w <- rep(1, nrow(x))
    ## covariance and centres as lists
    Cov <- list()
    for (g in levels(groups))
        Cov[[g]] <- cov.wt(x[groups == g,, drop = FALSE],
                           wt = w[groups == g])
    if (kind == "se")
        for(i in seq_len(length(Cov)))
            Cov[[i]]$cov <- Cov[[i]]$cov * sum(Cov[[i]]$wt^2)
    ## recycle colours
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = length(Cov))
    ## rgl::ellipse3d defaults to confidence envelopes, but we want to
    ## default to sd/se and only use confidence ellipses if conf is
    ## given
    if (missing(conf))
        t <- 1
    else
        t <- sqrt(qchisq(conf, 3))
    ## graph
    for(i in seq_len(length(Cov)))
        if (Cov[[i]]$n.obs > 3)
            plot3d(ellipse3d(Cov[[i]]$cov, centre = Cov[[i]]$center, t = t),
                   add = TRUE, alpha = alpha, col = col[i], ...)
}
