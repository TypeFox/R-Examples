##
## Copyright (c) 2009 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

panel.ablineq <-
    function(a = NULL, b = 0,
             h = NULL, v = NULL,
             reg = NULL, coef = NULL,
             pos = if (rotate) 1 else NULL,
             offset = 0.5, adj = NULL,
             at = 0.5,
             x = NULL, y = NULL,
             rotate = FALSE, srt = 0,
             label = NULL,
             varNames = alist(y = y, x = x),
             varStyle = "italic",
             fontfamily = "serif",
             digits = 3,
             r.squared = FALSE, sep = ", ", sep.end = "",
             col, col.text = add.text$col,
             col.line = add.line$col,
             ..., reference = FALSE)
{
    if (!is.null(label)) varStyle <- NULL
    ## work out colours for text and line
    add.text <- trellis.par.get("add.text")
    add.line <- if (reference)
        trellis.par.get("reference.line")
    else trellis.par.get("add.line")
    if (!missing(col) && missing(col.line))
        col.line <- col
    if (!missing(col) && missing(col.text))
        col.text <- col
    ## draw the line
    panel.abline(a = a, b = b, h = h, v = v, reg = reg, coef = coef,
                 col = col.line, ..., reference = reference)
    ## extract r.squared from model object if any
    if (!is.null(reg)) {
        a <- reg
    }
    if (isTRUE(r.squared)) {
        if (is.object(a) || is.list(a)) {
            r.squared <- round(summary(a)$r.squared, digits)
        } else {
            warning("r.squared = TRUE requires a model object")
        }
    }
    ## work out equation coefficients
    ## the following copied from lattice::panel.abline
    if (is.object(a) || is.list(a)) {
        p <- length(coefa <- as.vector(coef(a)))
        if (p > 2)
            warning("only using the first two of ", p, "regression coefficients")
        islm <- inherits(a, "lm")
        noInt <- if (islm)
            !as.logical(attr(stats::terms(a), "intercept"))
        else p == 1
        if (noInt) {
            a <- 0
            b <- coefa[1]
        }
        else {
            a <- coefa[1]
            b <- if (p >= 2)
                coefa[2]
            else 0
        }
    }
    if (!is.null(coef)) {
        if (!is.null(a))
            warning("'a' and 'b' are overridden by 'coef'")
        a <- coef[1]
        b <- coef[2]
    }
    if (length(h <- as.numeric(h)) > 0) {
        if (!is.null(a))
            warning("'a' and 'b' are overridden by 'h'")
        a <- h[1]
        b <- 0
    }
    if (length(a) > 1) {
        b <- a[2]
        a <- a[1]
    }
    ## construct the equation label
    if (length(as.numeric(v)) > 0) {
        ## vertical line (special case)
        if (!is.null(a))
            warning("'a' and 'b' are overridden by 'v'")
        if (is.null(x))
            x <- v[1]
        if (is.null(y))
            y <- convertY(unit(at, "npc"), "native", TRUE)
        if (is.null(label)) {
            v <- signif(v[1], digits)
            varNames <- c(as.list(varNames), v = v)
            label <- substitute(x == v, varNames)
        }
    } else {
        ## normal a+bx line
        if (is.null(x)) {
            ## work out start and end x values of visible line
            xlim <- current.panel.limits()$xlim
            ylim <- current.panel.limits()$ylim
            if (b == 0) {
                xx <- xlim
            } else {
                xx <- range((ylim - a) / b)
                xx <- pmin(max(xlim), xx)
                xx <- pmax(min(xlim), xx)
            }
            ## x position as fractional distance along line
            x <- min(xx) + at * abs(diff(xx))
        }
        if (is.null(y))
            y <- a + b * x
        if (is.null(label)) {
            a <- round(a, digits)
            b <- round(b, digits)
            varNames <- c(as.list(varNames), a = a, b = b)
            if (b == 0) {
                label <- substitute(y == a, varNames)
            } else if (a == 0) {
                label <- substitute(y == b * x, varNames)
            } else if (b > 0) {
                label <- substitute(y == a + b * x, varNames)
            } else {
                varNames$b <- abs(b)
                label <- substitute(y == a - b * x, varNames)
            }
        }
    }
    if (rotate) {
        if (length(as.numeric(v)) > 0) {
            srt <- 90
        } else {
            ## aspect ratio with respect to native coordinates
            asp <- with(lapply(current.panel.limits(), diff), ylim / xlim)
            ## aspect ratio of panel at *current* device size
            asp.cm <- with(lapply(current.panel.limits("cm"), diff), ylim / xlim)
            grad <- b * (asp.cm / asp)
            srt <- 180 * atan(grad) / pi
        }
    }
    if (is.numeric(r.squared)) {
        ## add R^2 = ... to label
        rsq.expr <- substitute(italic(R)^2 == z,
                            list(z = r.squared))
        label <- call("paste", label, sep, rsq.expr, sep.end)
    }
    ## wrap 'varStyle' function (or multiple functions) around label expression
    if (!is.null(varStyle)) {
        while (length(varStyle) > 0) {
            label <- call(varStyle[1], label)
            varStyle <- varStyle[-1]
        }
    }
    panel.text(x = x, y = y, labels = label,
               pos = pos, offset = offset, adj = adj,
               fontfamily = fontfamily, srt = srt,
               col = col.text, ...)
}

panel.lmlineq <-
    function(x, y, ...)
{
    if (length(x) > 1)
        panel.ablineq(lm(as.numeric(y) ~ as.numeric(x)), ...)
}

