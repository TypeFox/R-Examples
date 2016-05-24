blandAltman <- function(x, ...)
    UseMethod("blandAltman")

blandAltman.default <- function(x, y,
                                xlab = "Mean", ylab = "Difference",
                                main = "Bland-Altman plot",
                                sdLines = 2, ...) {
    xMean <- (x + y) / 2
    yDiff <- x - y
    
    ylim.min <- min(mean(yDiff, na.rm = TRUE) -
                        (sdLines + .1) * sd(yDiff, na.rm = TRUE),
                    min(yDiff, na.rm = TRUE), na.rm = TRUE)
    
    ylim.max <- max(mean(yDiff, na.rm = TRUE) +
                        (sdLines + 1) * sd(yDiff, na.rm = TRUE),
                    max(yDiff,na.rm = TRUE), na.rm = TRUE)
    
    plot(xMean, yDiff, xlab = xlab, ylab = ylab,
         main = main,
         ylim = c(ylim.min, ylim.max), ...)
    abline(h = mean(yDiff, na.rm = T) -
               c(-sdLines, 0, sdLines) * sd(yDiff, na.rm = T),
           lty = c(3, 2, 3))
}

blandAltman.formula <- function(formula, data, subset, na.action, ...) {
    # mostly copied from stats::t.test.formula
    if (missing(formula) ||
            (length(formula) != 3L) ||
            (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    
    m <- match.call(expand.dots = FALSE)
    
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    
    if (nlevels(g) != 2L) 
        stop("grouping factor must have exactly 2 levels")
    
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    do.call("blandAltman", c(DATA, list(...)))
}