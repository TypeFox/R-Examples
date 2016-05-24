diffPlot.formula <-
function (formula, data = NULL, plotFUN = mean, errFUN = c("ci", "se",
    "sd"), conf = 0.95, grp.names = NULL, var.equal = FALSE,
    paired = FALSE, ylim=NULL, ...)
{
    if (paired == TRUE)
        warning("You said paired=TRUE, but provided formula. Paired argument is ignored. For paired data, use 'x, y' notation.")
    se <- function(x) {
        x <- na.omit(x)
        res <- sqrt(var(x)/length(x))
        res
    }
    ci <- function(x) {
        x <- na.omit(x)
        alpha <- 1 - (1 - conf)/2
        res <- qt(alpha, length(x) - 2, lower.tail = T) * se(x)
        res
    }
    if (!is.null(errFUN)) {
        errFUN <- match.arg(errFUN)
    }
    dat <- model.frame(formula, data = data)
    dat <- na.omit(dat)
    if (length(unique(dat[, 2])) != 2) {
        stop("Grouping factor must have exactly two levels.")
    }
    res <- tapply(dat[, 1], dat[, 2], plotFUN)
    res <- c(res, res[2])
    Ns <- tapply(dat[, 1], dat[, 2], length)
    Vars <- tapply(dat[, 1], dat[, 2], var)
    if (var.equal == FALSE) {
        poolSE <- sqrt(sum(Vars/Ns))
        df <- (sum(Vars/Ns)^2)/(sum(Vars^2/(Ns^2 * (Ns - 1))))
    }
    if (var.equal == TRUE) {
        poolSE <- sqrt(sum((Ns - 1) * Vars)/(sum(Ns - 1))) *
            sqrt(sum(1/Ns))
        df <- sum(Ns - 1)
    }
    diffBar <- ifelse(errFUN == "ci", qt(1 - (1 - conf)/2, df,
        lower.tail = T) * poolSE, ifelse(errFUN == "se", poolSE,
        ifelse(errFUN == "sd", sqrt(sum(Vars)), 0)))
    e <- c(tapply(dat[, 1], dat[, 2], errFUN), diffBar)
    e <- ifelse(res <= 0, -e, e)
    if(!is.null(ylim)) {lims <- ylim}
    else {
      lims <- c(min(res + e, res - e) - abs(0.4 * min(res + e,
          res - e)), max(res + e, res - e) + 0.4 * abs(max(res +
          e, res - e)))
    }
    par(mar = c(5, 6, 4, 5), las = 1)
    plot(c(1, 2, 4), res, pch = c(19, 19, 17), xaxt = "n", xlim = c(0.4,
        0.4 + 4), ylim = c(lims), bty = "l", ...)
    if (!is.null(errFUN)) {
        arrows(c(1, 2, 4), res + e, c(1, 2, 4), res - e, angle = 90,
            code = 3, length = 0.08)
        arrows(2, res[2], 4.4, res[3], angle = 90, code = 3,
            length = 0, lty = 2)
        arrows(1, res[1], 4.4, res[1], angle = 90, code = 3,
            length = 0, lty = 2)
        if (res[1] >= res[2]) {
            val <- axTicks(4) - axTicks(4)[min(which(axTicks(4) >
                median(axTicks(4))))]
            loc <- axTicks(4) - (axTicks(4)[min(which(axTicks(4) >
                median(axTicks(4))))] - res[1])
        }
        if (res[1] < res[2]) {
            val <- axTicks(4) - axTicks(4)[max(which(axTicks(4) <
                median(axTicks(4))))]
            loc <- axTicks(4) - (axTicks(4)[max(which(axTicks(4) <
                median(axTicks(4))))] - res[1])
        }
        axis(4, at = loc, labels = val, ...)
    }
    if (is.null(grp.names)) {
        grp.names <- unique(dat[, 2])
    }
    axis(1, at = c(1, 2, 4), labels = c(grp.names, "Difference"))
}
