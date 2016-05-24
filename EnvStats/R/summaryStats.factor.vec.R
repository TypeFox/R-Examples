summaryStats.factor.vec <-
function (x, combine.levels = TRUE, digits = max(3, getOption("digits") - 
    3), digit.type = "round", show.na = TRUE, show.0.na = FALSE) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    levels.x <- levels(x)
    n.levels <- length(levels.x)
    n <- summary(x)
    names.n <- names(n)
    if (length(n) == n.levels) {
        Pct <- do.call(digit.type, args = list(x = 100 * n/sum(n), 
            digits = digits))
        if (combine.levels) {
            n <- c(n, sum(n))
            Pct <- c(Pct, 100)
            names.n <- c(names.n, "Combined")
        }
        if (show.na && show.0.na) {
            n <- c(n, 0)
            Pct <- c(Pct, NA)
            names.n <- c(names.n, "NA's")
        }
        ret.n <- n
    }
    else {
        ret.n <- n[levels.x]
        names.n <- names(ret.n)
        Pct <- do.call(digit.type, args = list(x = 100 * ret.n/sum(ret.n), 
            digits = digits))
        if (combine.levels) {
            ret.n <- c(ret.n, sum(ret.n))
            Pct <- c(Pct, 100)
            names.n <- c(names.n, "Combined")
        }
        if (show.na) {
            ret.n <- c(ret.n, n["NA's"])
            Pct <- c(Pct, NA)
            names.n <- c(names.n, "NA's")
        }
    }
    mat <- cbind(N = ret.n, Pct)
    dimnames(mat)[[1]] <- names.n
    mat
}
