DanielPlot <-
function (fit, code = FALSE, faclab = NULL, block = FALSE, datax = TRUE, 
    half = FALSE, pch = "*", cex.fac = par("cex.lab"), cex.lab = par("cex.lab"), 
    cex.pch = par("cex.axis"), ...) 
{
    if (any(names(coef(fit)) == "(Intercept)")) {
        factor.effects <- 2 * coef(fit)[-1]
    }
    else {
        factor.effects <- 2 * coef(fit)
    }
    names(factor.effects) <- attr(fit$terms, "term.labels")
    factor.effects <- factor.effects[!is.na(factor.effects)]
    if (half) {
        tn <- data.frame(x = qnorm(0.5 * ((rank(abs(factor.effects)) - 
            0.5)/length(factor.effects) + 1)), x = abs(factor.effects))
        names(tn$x) <- names(factor.effects)
        xlab <- "half-normal score"
        ylab <- "absolute effects"
    }
    else {
        tn <- qqnorm(factor.effects, plot = FALSE)
        xlab <- "normal score"
        ylab <- "effects"
    }
    if (datax) {
        tmp <- tn$x
        tn$x <- tn$y
        tn$y <- tmp
        tmp <- xlab
        xlab <- ylab
        ylab <- tmp
    }
    labx <- names(factor.effects)
    laby <- 1:length(tn$y)
    points.labels <- names(factor.effects)
    plot.default(tn, xlim = c(min(tn$x), max(tn$x) + diff(range(tn$x))/5), 
        pch = pch, xlab = xlab, ylab = ylab, cex.lab = cex.lab, 
        ...)
    if (is.null(faclab)) {
        if (!code) {
            effect.code <- labx
        }
        else {
            terms.ord <- attr(fit$terms, "order")
            max.order <- max(terms.ord)
            no.factors <- length(terms.ord[terms.ord == 1])
            factor.label <- attr(fit$terms, "term.labels")[1:no.factors]
            factor.code <- LETTERS[1:no.factors]
            if (block) 
                factor.code <- c("BK", factor.code)
            texto <- paste(factor.code[1], "=", factor.label[1])
            for (i in 2:no.factors) {
                texto <- paste(texto, ", ", factor.code[i], "=", 
                  factor.label[i])
            }
            mtext(side = 1, line = 2.5, texto, cex = cex.fac)
            get.sep <- function(string, max.order) {
                k <- max.order - 1
                get.sep <- rep(0, k)
                j <- 1
                for (i in 1:nchar(string)) {
                  if (substring(string, i, i) == ":") {
                    get.sep[j] <- i
                    if (j == k) 
                      break
                    j <- j + 1
                  }
                }
                get.sep
            }
            labeling <- function(string, get.sep, max.order, 
                factor.code, factor.label) {
                labeling <- ""
                sep <- get.sep(string, max.order)
                sep <- sep[sep > 0]
                n <- length(sep) + 1
                if (n > 1) {
                  sep <- c(0, sep, nchar(string) + 1)
                  for (i in 1:n) {
                    labeling <- paste(labeling, sep = "", factor.code[factor.label == 
                      substring(string, sep[i] + 1, sep[i + 1] - 
                        1)][1])
                  }
                }
                else labeling <- paste(labeling, sep = "", factor.code[factor.label == 
                  string][1])
                labeling
            }
            effect.code <- rep("", length(terms.ord))
            for (i in 1:length(terms.ord)) {
                effect.code[i] <- labeling(names(tn$x)[i], get.sep, 
                  max.order, factor.code, factor.label)
            }
        }
        text(tn, paste("   ", effect.code), cex = cex.pch, adj = 0, 
            xpd = NA)
    }
    else {
        if (!is.list(faclab)) 
            stop("* Argument 'faclab' has to be NULL or a list with idx and lab objects")
        text(tn$x[faclab$idx], tn$y[faclab$idx], labels = faclab$lab, 
            cex = cex.fac, adj = 0)
    }
    invisible(cbind(as.data.frame(tn), no = 1:length(tn$x)))
}
