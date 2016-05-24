plotCor <- function(tab,
                    method = c("pearson", "kendall", "spearman"),
                    formatCor, colCor, ...) {
    if (missing(formatCor)) {
        formatCor <- function(test, n = 3, ...) {
            if (is.null(test$conf.int)) {
                paste(round(test$estimate, n),
                      if (test$p.value < .001) "p < 0.001" else sprintf("p = %.03f", test$p.value),
                      sep = "\n")
            } else {
                paste(round(test$estimate, n),
                      sprintf(sprintf("[%%.%if ; %%.%if]", n, n),
                              test$conf.int[1],
                              test$conf.int[2]),
                      if (test$p.value < .001) "p < 0.001" else sprintf("p = %.03f", test$p.value),
                      sep = "\n")
            }
        }
    }
    
    if (missing(colCor)) {
        colCor <- function(test) {
            rgb(red = (r <- test$estimate > 0),
                green = 0,
                blue = ! r,
                alpha = abs(test$estimate))
        }
    }
    
    nc <- ncol(tab)
    stopifnot(nc > 1)
    method <- match.arg(method)
    op <- par(mfrow = c(nc, nc),
              mai = c(0, 0, 0, 0),
              xaxt = "n",
              yaxt = "n")
    
    try({
        for (i in 1:nc) {
            for (j in 1:nc) {
                if (i == j) {
                    plot.new()
                    text(.5, .5, names(tab)[i])
                } else {
                    if (j > i) {
                        test <- cor.test(tab[, j],
                                         tab[, i],
                                         method = method)
                        test.text <- formatCor(test, ...)
                        plot.new()
                        rect(0, 0, 1, 1, col = colCor(test),
                             border = NA)
                        text(.5, .5, test.text)
                    } else {
                        plot(tab[, j], tab[, i])
                    }
                }
            }
        }})
    
    par(op)
}

plotDensity <- function(x, col = "black", lwd = 2, ...) {
    d <- density(na.omit(x))
    h <- hist(x, plot = FALSE)
    plot(h, xlim = c(min(d$x), max(d$x)),
         ylim = c(0, max(c(d$y, h$density))),
         freq = FALSE, ...)
    lines(d, lwd = lwd, col = col)
}

plotICC <- function(x, ...)
    UseMethod("plotICC")

plotICC.default <- function(x, subjects, p = FALSE, ptype = 4, ...) {
    lMax <- max(tabS <- table(subjects)) + 1
    id <- names(tabS)
    FunC <- function(a, lMax) {
        c(a, rep(NA, lMax - length(a)))
    }
    listY <- lapply(tapply(x, subjects, c), FunC, lMax = lMax)
    matY <- matrix(unlist(listY), nrow = lMax)
    matY <- matY[ , order(colMeans(matY, na.rm = TRUE))]
    matX <- matrix(rep(seq_along(id), each = lMax), nrow = lMax)
    plot(matX, matY, type = "l", ...)
    if (p) points(matX, matY, pch = ptype)
}

plotICC.data.frame <- function(x, p = FALSE, ptype = 4, ...) {
    nObs <- ncol(x)
    nSubjects <- nrow(x)
    x <- t(x)
    dim(x) <- NULL
    subjects <- rep(seq(length.out = nSubjects), each = nObs)
    plotICC(x, subjects, p = p, ptype = ptype, ...)
}

plotICC.matrix <- function(x, p = FALSE, ptype = 4, ...) {
    x <- as.data.frame(x)
    plotICC(x, p = p, ptype = ptype, ...)
}

plotScatter <- function(x, y, lty = 1, lwd = 2, colLine = "red", ...) {
    m <- na.omit(cbind(x, y))
    plot(m, ...)
    lines(lowess(m), col = colLine, lwd = lwd, lty = lty)
}

qqplot2 <- function(y, fQuant = function(q, x) qnorm(q,
                                                     mean(x, na.rm = TRUE),
                                                     sd(x, na.rm = TRUE)), 
                    line = TRUE, xlab = "Theoretical Quantiles", ...) {
    qy <- ppoints(length(y))
    qqplot(fQuant(qy, y), y,
           xlab = xlab, ...)
    if (line)
        abline(a = 0, b = 1, lty = 2)
}

ablineCI <- function(x, level = .95, lty = 2, ...) {
    p <- length(cx <- coef(x))
    if (p > 2) 
        warning(sprintf("only using the first two of %d regression coefficients", 
                        p))
    stopifnot(inherits(x, "lm"))
    
    ci <- predict(x, interval = "confidence", level = level)
    
    xname <- (n <- names(cx))[n != "(Intercept)"][1]
    xval <- x$model[[xname]]
    
    ci <- ci[order(xval), ]
    xval <- sort(xval)
    
    lines(xval, ci[, "lwr"], lty = lty, ...)
    lines(xval, ci[, "upr"], lty = lty, ...)
    
    invisible(structure(cbind(xval, ci),
                        dimnames = list(NULL, c(xname, colnames(ci)))))
}
