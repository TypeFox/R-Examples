ccdplot <-
function (x, remove.absolute = NA, remove.ratio = NA, drawcomposite = TRUE, 
    jump = NA, xlab = "Observations", ylab = "Cumulatives", ...) 
{
    check.remove.params(remove.absolute, remove.ratio)
    if (!is.list(x)) {
        x.removed <- remove.vector(x, remove.absolute, remove.ratio)
        x.sortsum <- calculate.sortsum(x.removed)
        plot(1:length(x.sortsum), x.sortsum, type = "l", xlab = xlab, 
            ylab = ylab, ...)
        abline(h = 0, lty = "dotted")
    }
    else {
        x.removed <- remove.list(x, remove.absolute, remove.ratio)
        characteristics.number <- length(x.removed)
        x.composite <- c()
        x.sortsum <- list()
        absoluteMin <- Inf
        absoluteMax <- -Inf
        for (i in 1:characteristics.number) {
            x.sortsum[[i]] <- calculate.sortsum(x.removed[[i]])
            x.composite <- c(x.composite, x.removed[[i]])
            actualMin <- min(x.sortsum[[i]])
            if (!is.na(actualMin) && actualMin < absoluteMin) {
                absoluteMin <- actualMin
            }
            actualMax <- max(x.sortsum[[i]])
            if (!is.na(actualMax) && actualMax > absoluteMax) {
                absoluteMax <- actualMax
            }
        }
        x.composite.sortsum <- calculate.sortsum(x.composite)
        actualMin <- min(x.composite.sortsum)
        if (actualMin < absoluteMin) {
            absoluteMin <- actualMin
        }
        actualMax <- max(x.composite.sortsum)
        if (actualMax > absoluteMax) {
            absoluteMax <- actualMax
        }
        if (is.na(jump)) {
            jump <- length(x.composite.sortsum)/50
        }
        xmax <- length(x.composite.sortsum)
        if (!drawcomposite) {
            xmax <- xmax + (characteristics.number - 1) * jump
        }
        plot(c(1, xmax), c(absoluteMin, absoluteMax), xlab = xlab, 
            ylab = ylab, type = "n", ...)
        abline(h = 0, lty = "dotted")
        cumulativeLength <- 0
        loopend <- characteristics.number
        if (drawcomposite) {
            points(1:length(x.composite.sortsum), x.composite.sortsum, 
                type = "l", ...)
            loopend <- loopend - 1
        }
        for (i in 1:loopend) {
            x.from <- i * jump + 1 + cumulativeLength
            x.to <- i * jump + cumulativeLength + length(x.sortsum[[i]])
            if (!drawcomposite) {
                x.from <- x.from - jump
                x.to <- x.to - jump
            }
            points(x.from:x.to, x.sortsum[[i]], type = "l", ...)
            cumulativeLength <- cumulativeLength + length(x.sortsum[[i]])
        }
        if (drawcomposite) {
            points((-jump + 1 + cumulativeLength):(-jump + cumulativeLength + 
                length(x.sortsum[[characteristics.number]])), 
                x.sortsum[[characteristics.number]], type = "l", 
                ...)
        }
    }
}
