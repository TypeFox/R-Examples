# $Id: qqnorm.aov.R 625 2005-06-09 14:20:30Z nj7w $

if(is.R())
qqnorm.aov <- function (y, full = FALSE, label = FALSE, omit = NULL,
                        xlab = paste(if(full) "" else "Half", " Normal plot"),
                        ylab = "Effects", ...)
{
    r <- y$rank
    eff <- if (full)
        effects(y, set.sign = TRUE)[1:r]
    else abs(effects(y))[1:r]
    na <- names(eff)
    int <- match("(Intercept)", na)
    if (!is.null(omit)) {
        if (is.character(omit)) {
            int <- c(int, match(omit, na))
        }
        else int <- c(int, omit)
    }
    int <- int[!is.na(int)]
    if (length(int))
        eff <- eff[-int]
    n <- length(eff)
    if (n <= 0)
        stop("Not enough effects")
    ord <- order(eff)
    na <- names(eff)
    P <- if (full)
        ppoints(n)
    else ((1:n) + n)/(2 * n + 1)
    Q <- qnorm(P)
    plot(x = Q, y = eff[ord], xlab = xlab, ylab = ylab, ...)
    if (label && dev.interactive())
        identify(Q, eff[ord], names(eff)[ord])
}


