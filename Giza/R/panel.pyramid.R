panel.pyramid <-
function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio), 
    horizontal = TRUE, origin = NULL, reference = TRUE, stack = TRUE, 
    groups = NULL, col = if (is.null(groups)) plot.polygon$col else superpose.polygon$col, 
    border = if (is.null(groups)) plot.polygon$border else superpose.polygon$border, 
    lty = if (is.null(groups)) plot.polygon$lty else superpose.polygon$lty, 
    lwd = if (is.null(groups)) plot.polygon$lwd else superpose.polygon$lwd, 
    ..., identifier = "barchart") 
{
    plot.polygon <- trellis.par.get("plot.polygon")
    superpose.polygon <- trellis.par.get("superpose.polygon")
    reference.line <- trellis.par.get("reference.line")
    keep <- (function(x, y, groups, subscripts, ...) {
        !is.na(x) & !is.na(y) & if (is.null(groups)) 
            TRUE
        else !is.na(groups[subscripts])
    })(x = x, y = y, groups = groups, ...)
    if (!any(keep)) 
        return()
    x <- as.numeric(x[keep])
    y <- as.numeric(y[keep])
    if (!is.null(groups)) {
        groupSub <- function(groups, subscripts, ...) groups[subscripts[keep]]
        if (!is.factor(groups)) 
            groups <- factor(groups)
        nvals <- nlevels(groups)
        groups <- as.numeric(groupSub(groups, ...))
    }
    if (horizontal) {
        if (is.null(groups)) {
            if (is.null(origin)) {
                origin <- current.panel.limits()$xlim[1]
                reference <- FALSE
            }
            height <- box.width
            if (reference) 
                panel.abline(v = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            panel.rect(x = rep(origin, length(y)), y = y, height = rep(height, 
                length(y)), width = x - origin, border = border, 
                col = col, lty = lty, lwd = lwd, just = c("left", 
                  "centre"), identifier = identifier)
        }
        else if (stack) {
            if (!is.null(origin) && origin != 0) 
                warning("'origin' forced to 0 for stacked bars")
            col <- rep(col, length.out = nvals)
            border <- rep(border, length.out = nvals)
            lty <- rep(lty, length.out = nvals)
            lwd <- rep(lwd, length.out = nvals)
            height <- box.width
            if (reference) 
                panel.abline(v = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            for (i in unique(y)) {
                ok <- y == i
                ord <- sort.list(groups[ok])
                pos <- x[ok][ord] > 0
                nok <- sum(pos, na.rm = TRUE)
                if (nok > 0) 
                  panel.rect(x = cumsum(c(0, x[ok][ord][pos][-nok])), 
                    y = rep(i, nok), col = col[groups[ok][ord][pos]], 
                    border = border[groups[ok][ord][pos]], lty = lty[groups[ok][ord][pos]], 
                    lwd = lwd[groups[ok][ord][pos]], height = rep(height, 
                      nok), width = x[ok][ord][pos], just = c("left", 
                      "centre"), identifier = paste(identifier, 
                      "pos", i, sep = "."))
                neg <- x[ok][ord] < 0
                nok <- sum(neg, na.rm = TRUE)
                if (nok > 0) 
                  panel.rect(x = cumsum(c(0, x[ok][ord][neg][-nok])), 
                    y = rep(i, nok), col = col[groups[ok][ord][neg]], 
                    border = border[groups[ok][ord][neg]], lty = lty[groups[ok][ord][neg]], 
                    lwd = lwd[groups[ok][ord][neg]], height = rep(height, 
                      nok), width = x[ok][ord][neg], just = c("left", 
                      "centre"), identifier = paste(identifier, 
                      "neg", i, sep = "."))
            }
        }
        else {
            if (is.null(origin)) {
                origin <- current.panel.limits()$xlim[1]
                reference <- FALSE
            }
            col <- rep(col, length.out = nvals)
            border <- rep(border, length.out = nvals)
            lty <- rep(lty, length.out = nvals)
            lwd <- rep(lwd, length.out = nvals)
            height <- box.width/nvals
            if (reference) 
                panel.abline(v = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            for (i in unique(y)) {
                ok <- y == i
                nok <- sum(ok, na.rm = TRUE)
                panel.rect(x = rep(origin, nok), y = (i + height * 
                  (groups[ok] - (nvals + 1)/2)), col = col[groups[ok]], 
                  border = border[groups[ok]], lty = lty[groups[ok]], 
                  lwd = lwd[groups[ok]], height = rep(height, 
                    nok), width = x[ok] - origin, just = c("left", 
                    "centre"), identifier = identifier)
            }
        }
    }
    else {
        if (is.null(groups)) {
            if (is.null(origin)) {
                origin <- current.panel.limits()$ylim[1]
                reference <- FALSE
            }
            width <- box.width
            if (reference) 
                panel.abline(h = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            panel.rect(x = x, y = rep(origin, length(x)), col = col, 
                border = border, lty = lty, lwd = lwd, width = rep(width, 
                  length(x)), height = y - origin, just = c("centre", 
                  "bottom"), identifier = identifier)
        }
        else if (stack) {
            if (!is.null(origin) && origin != 0) 
                warning("'origin' forced to 0 for stacked bars")
            col <- rep(col, length.out = nvals)
            border <- rep(border, length.out = nvals)
            lty <- rep(lty, length.out = nvals)
            lwd <- rep(lwd, length.out = nvals)
            width <- box.width
            if (reference) 
                panel.abline(h = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            for (i in unique(x)) {
                ok <- x == i
                ord <- sort.list(groups[ok])
                pos <- y[ok][ord] > 0
                nok <- sum(pos, na.rm = TRUE)
                if (nok > 0) 
                  panel.rect(x = rep(i, nok), y = cumsum(c(0, 
                    y[ok][ord][pos][-nok])), col = col[groups[ok][ord][pos]], 
                    border = border[groups[ok][ord][pos]], lty = lty[groups[ok][ord][pos]], 
                    lwd = lwd[groups[ok][ord][pos]], width = rep(width, 
                      nok), height = y[ok][ord][pos], just = c("centre", 
                      "bottom"), identifier = paste(identifier, 
                      "pos", i, sep = "."))
                neg <- y[ok][ord] < 0
                nok <- sum(neg, na.rm = TRUE)
                if (nok > 0) 
                  panel.rect(x = rep(i, nok), y = cumsum(c(0, 
                    y[ok][ord][neg][-nok])), col = col[groups[ok][ord][neg]], 
                    border = border[groups[ok][ord][neg]], lty = lty[groups[ok][ord][neg]], 
                    lwd = lwd[groups[ok][ord][neg]], width = rep(width, 
                      nok), height = y[ok][ord][neg], just = c("centre", 
                      "bottom"), identifier = paste(identifier, 
                      "neg", i, sep = "."))
            }
        }
        else {
            if (is.null(origin)) {
                origin <- current.panel.limits()$ylim[1]
                reference = FALSE
            }
            col <- rep(col, length.out = nvals)
            border <- rep(border, length.out = nvals)
            lty <- rep(lty, length.out = nvals)
            lwd <- rep(lwd, length.out = nvals)
            width <- box.width/nvals
            if (reference) 
                panel.abline(h = origin, col = reference.line$col, 
                  lty = reference.line$lty, lwd = reference.line$lwd, 
                  identifier = paste(identifier, "abline", sep = "."))
            for (i in unique(x)) {
                ok <- x == i
                nok <- sum(ok, na.rm = TRUE)
                panel.rect(x = (i + width * (groups[ok] - (nvals + 
                  1)/2)), y = rep(origin, nok), col = col[groups[ok]], 
                  border = border[groups[ok]], lty = lty[groups[ok]], 
                  lwd = lwd[groups[ok]], width = rep(width, nok), 
                  height = y[ok] - origin, just = c("centre", 
                    "bottom"), identifier = identifier)
            }
        }
    }
}
