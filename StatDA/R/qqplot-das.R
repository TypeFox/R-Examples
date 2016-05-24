"qqplot.das" <-
function (x, distribution = "norm", ylab = deparse(substitute(x)), 
    xlab = paste(distribution, "quantiles"), main = "", las = par("las"), 
    datax=FALSE, envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, 
    pch = 1, line = c("quartiles", "robust", "none"), cex=1,xaxt="s",
    add.plot=FALSE, xlim = NULL, ylim=NULL, ...) 
{
# Modified by PF, 06.11.2009
# original function: qq.plot from library(car)
# Modification: additional argument "datax"
# If datax=TRUE: x and y axes are exchanged
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    d.function <- eval(parse(text = paste("d", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, ...)
    if(!add.plot) {
      if (datax)
        plot(ord.x, z, xlab = ylab, ylab = xlab, main = main,
        las = las, col = col, pch = pch, cex = cex, xaxt = xaxt,
        ylim = ylim,xlim=xlim)
      else plot(z, ord.x, xlab = xlab, ylab = ylab, main = main,
          las = las, col = col, pch = pch, cex = cex, xaxt = xaxt, 
          ylim = ylim,xlim=xlim)
    }
    else {
      if (datax)
        points(ord.x, z, xlab = ylab, ylab = xlab, main = main,
        las = las, col = col, pch = pch, cex = cex, xaxt = xaxt)
      else points(z, ord.x, xlab = xlab, ylab = ylab, main = main,
        las = las, col = col, pch = pch, cex = cex, xaxt = xaxt)

    }

    if (line == "quartiles") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75), ...)
        if (datax) {
          b <- (Q.z[2] - Q.z[1])/(Q.x[2] - Q.x[1])
          a <- Q.z[1] - b * Q.x[1]
        }
        else {
          b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
          a <- Q.x[1] - b * Q.z[1]
        }
        abline(a, b, col = col, lwd = lwd)
    }
    if (line == "robust") {
        #if (!require("MASS")) 
        #    stop("MASS package not available")
        if (datax)
          coef <- coefficients(rlm(z ~ ord.x))
        else coef <- coefficients(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
        abline(a, b)
    }
    if (line != "none" & envelope != FALSE) {
        zz <- qnorm(1 - (1 - envelope)/2)
        if (datax) {
          b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
          a <- Q.x[1] - b * Q.z[1]
          SE <- (b/d.function(z, ...)) * sqrt(P * (1 - P)/n)
          fit.value <- a + b * z
          upper <- fit.value + zz * SE
          lower <- fit.value - zz * SE
          lines(upper, z, lty = 2, lwd = lwd/2, col = col)
          lines(lower, z, lty = 2, lwd = lwd/2, col = col)
        }
        else {
          SE <- (b/d.function(z, ...)) * sqrt(P * (1 - P)/n)
          fit.value <- a + b * z
          upper <- fit.value + zz * SE
          lower <- fit.value - zz * SE
          lines(z, upper, lty = 2, lwd = lwd/2, col = col)
          lines(z, lower, lty = 2, lwd = lwd/2, col = col)
        }
    }
    if (labels[1] == TRUE & length(labels) == 1) 
        labels <- seq(along = z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along = x)[good][ord][selected]
    }
    if (is.null(result)) 
        invisible(result)
    else sort(result)
}
