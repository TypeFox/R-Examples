biserial.cor <-
function (x, y, use = c("all.obs", "complete.obs"), level = 1) {
    if (!is.numeric(x))
        stop("'x' must be a numeric variable.\n")
    y <- as.factor(y)
    if (length(levs <- levels(y)) > 2)
        stop("'y' must be a dichotomous variable.\n")
    if (length(x) != length(y))
        stop("'x' and 'y' do not have the same length")
    use <- match.arg(use)
    if (use == "complete.obs") {
        cc.ind <- complete.cases(x, y)
        x <- x[cc.ind]
        y <- y[cc.ind]
    } 
    ind <- y == levs[level]
    diff.mu <- mean(x[ind]) - mean(x[!ind])
    prob <- mean(ind)
    diff.mu * sqrt(prob * (1 - prob)) / sd(x)
}
