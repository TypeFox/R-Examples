`deshrinkPred` <- function(x, coef,
                           type = c("inverse", "classical",
                           "expanded", "none","monotonic")) {
    if(missing(type))
        type <- "inverse"
    type <- match.arg(type)
    res <- switch(type,
                  inverse = coef[1] + (coef[2] * x),
                  classical = (x - coef[1]) / coef[2],
                  expanded = coef[1] + (coef[2] * x),
                  none = coef[1] + (coef[2] * x),
                  monotonic = drop(Predict.matrix(coef$sm,
                  data.frame(wa.env = x)) %*% coef$p))
    return(res)
}
