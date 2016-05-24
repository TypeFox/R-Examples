
#' Augmented functions
#' 
#' These functions all behave similarly to the functions with the initial
#' \code{x} removed from their names.
#' 
#' 
#' @aliases xplot xplot.default xplot.lm
#' @param \dots arguments passed to other functions.
#' @param x A numeric vector or formula.
#' @param which A numeric vector indicating which plots to produce
#' @param panel.default default panel function
#' @param sub.caption secondary caption
#' @param main as in \code{xyplot}
#' @param type as in \code{xyplot}
#' @param pch as in \code{xyplot}
#' @param lty as in \code{xyplot}
#' @param qqline a logical
#' @param caption caption for the plot
#' @param print.plots a logical
#' @param ask a logical
#' @param addline.col color for added lines
#' @param line.col color for lines
#' @param symbol.col color for symbols
#' @param id.n a numeric
#' @param labels.id a character vector of labels 
#' @param cex.id cex for ids
#' @param cook.levels a logical
#' @param add.smooth a logical
#' @param label.pos position for labels, one of \code{"left"} or \code{"right"}
#' @param cex.caption cex for the caption
#' @param panel a panel function
#' @seealso \code{\link{plot}}.
#' @examples
#' 
#' x <- runif(20)
#' xplot( lm ( 2*x + 5 + rnorm(20) ~ x ) )
#' 

#' @export
xplot <-
function (x, ...) 
{
    UseMethod("xplot")
}

#' @rdname xplot
#' @method xplot default
#' @export
xplot.default <-
function (...) 
{
    plot(...)
}


#' @rdname xplot
#' @method xplot lm
#' @export
xplot.lm <-
function (x, which = c(1L:3, 5), 
		  caption = captions,
	panel.default = if (add.smooth) panel.xyplotsmooth else panel.xyplotpoints, 
    sub.caption = NULL, main = "", print.plots = TRUE, ask = 1 < 
        length(which) && dev.interactive(), type = "p", pch = trellis.par.get("plot.symbol")$pch, 
    addline.col = trellis.par.get("add.line")$col, line.col = trellis.par.get("plot.line")$col, 
    symbol.col = trellis.par.get("plot.symbol")$col, lty = trellis.par.get("superpose.line")$lty, 
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.7, 
    qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = TRUE, 
    label.pos = c("left", "right"), cex.caption = 1) 
{
	# if you don't like this style, complain to CRAN about 90-character line limits or to
	# roxygen about inserting natural line breaks.
	captions <- list("Residuals vs Fitted", 
    				 "Normal Q-Q", 
					 "Scale-Location", 
					 "Cook's distance", 
					 "Residuals vs Leverage", 
                     "Cook's distance vs Leverage")
    lty = rep(lty, 5)
    old.theme <- trellis.par.get()
    trellis.par.set(list(superpose.line = list(lty = lty), add.line = list(col = addline.col), 
        plot.line = list(col = line.col), plot.symbol = list(col = symbol.col), 
        par.main.text = list(cex = 0.9), par.sub.text = list(cex = 0.8)))
    results <- list()
    dropInf <- function(x, h) {
        if (missing(h)) {
            return(x)
        }
        if (any(isInf <- h >= 1)) {
            warning("Not plotting observations with leverage one:\n  ", 
                paste(which(isInf), collapse = ", "), call. = FALSE)
            x[isInf] <- NaN
        }
        x
    }
    if (!inherits(x, "lm")) 
        stop("use only with \"lm\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
        stop("'which' must be in 1L:6")
    isGlm <- inherits(x, "glm")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    r <- residuals(x)
    yh <- predict(x)
    w <- weights(x)
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    n <- length(r)
    if (any(show[2L:6L])) {
        s <- if (inherits(x, "rlm")) 
            x$s
        else if (isGlm) 
            sqrt(summary(x)$dispersion)
        else sqrt(deviance(x)/df.residual(x))
        hii <- lm.influence(x, do.coef = FALSE)$hat
        if (any(show[4L:6L])) {
            cook <- if (isGlm) 
                cooks.distance(x)
            else cooks.distance(x, sd = s, res = r)
        }
    }
    if (any(show[2L:3L])) {
        ylab23 <- if (isGlm) 
            "Std. deviance resid."
        else "Standardized residuals"
        r.w <- if (is.null(w)) 
            r
        else sqrt(w) * r
        rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    }
    if (any(show[5L:6L])) {
        r.hat <- range(hii, na.rm = TRUE)
        isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 * 
            mean(hii, na.rm = TRUE)
    }
    if (any(show[c(1L, 3L)])) 
        l.fit <- if (isGlm) 
            "Predicted values"
        else "Fitted values"
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0L || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), 
                domain = NA)
    }
    if (id.n > 0L) {
        if (is.null(labels.id)) 
            labels.id <- paste(1L:n)
        iid <- 1L:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        if (any(show[2L:3L])) 
            show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    }
    getCaption <- function(k) as.graphicsAnnot(unlist(caption[k]))
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2L] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1L], "c")
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr) 
            paste(substr(cc[1L], 1L, min(75L, nc)), "...")
        else cc[1L]
    }
    one.fig <- 1 < length(which)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (show[1L]) {
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0) 
            ylim <- extendrange(r = ylim, f = 0.08)
        newplot <- xyplot(r ~ yh, xlab = l.fit, ylab = "Residuals", 
            ylim = ylim, type = "n", main = if (one.fig) 
                getCaption(1)
            else NULL, panel = function(x, y) {
                panel.abline(h = 0, lwd = 2, lty = trellis.par.get("add.line")$lty, 
                  col = trellis.par.get("add.line")$col)
                panel.xyplot(x, y, ...)
                if (id.n > 0) {
                  grid.identify.points(x, y, show.r, cex = cex.id)
                }
            }, ...)
        results <- c(results, list(newplot))
    }
    if (show[2L]) {
        ylim <- range(rs, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        newplot <- qqmath(rs, ylab = ylab23, main = if (one.fig) 
            getCaption(2)
        else NULL, sub = if (one.fig) 
            sub.caption
        else NULL, panel = function(x, ...) {
            if (qqline) {
                panel.qqmathline(x, lty = trellis.par.get("add.line")$lty, 
                  lwd = 2, ...)
            }
            panel.qqmath(x, pch = 16, ...)
            if (id.n > 0) 
                grid.identify.points(x = qnorm(ppoints(length(x))[rank(x)]), 
                  y = x, show.rs, cex = cex.id)
        }, ...)
        results <- c(results, list(newplot))
    }
    if (show[3L]) {
        sqrtabsr <- sqrt(abs(rs))
        ylim <- c(0, 1.05 * max(sqrtabsr, na.rm = TRUE))
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
        yhn0 <- if (is.null(w)) 
            yh
        else yh[w != 0]
        newplot <- xyplot(sqrtabsr ~ yhn0, xlab = l.fit, ylab = yl, 
            main = if (one.fig) 
                getCaption(3)
            else NULL, sub = if (one.fig) 
                sub.caption
            else NULL, ylim = ylim, panel = function(x, y, ...) {
                panel.default(x, y, line.col = line.col, ...)
                if (id.n > 0) {
                  grid.identify.points(yhn0, sqrtabsr, show.rs, 
                    cex = cex.id)
                }
            }, ...)
        results <- c(results, list(newplot))
    }
    if (show[4L]) {
        if (id.n > 0) {
            show.r <- order(-cook)[iid]
            ymx <- cook[show.r[1L]] * 1.075
        }
        else ymx <- max(cook, na.rm = TRUE)
        newplot <- xyplot(cook ~ 1:length(cook), type = "h", 
            xlab = "Obs. number", ylab = "Cook's distance", main = if (one.fig) 
                getCaption(4)
            else NULL, sub = if (one.fig) 
                sub.caption
            else NULL, panel = function(x, y, ...) {
                panel.default(x, y, line.col = line.col, ...)
                if (id.n > 0) {
                  grid.identify.points(x, y, show.r, adj.x = FALSE, 
                    adj.y = T, cex = cex.id)
                }
            }, ...)
        results <- c(results, list(newplot))
    }
    if (show[5L]) {
        ylab5 <- if (isGlm) 
            "Std. Pearson resid."
        else "Standardized residuals"
        r.w <- residuals(x, "pearson")
        if (!is.null(w)) 
            r.w <- r.w[wind]
        rsp <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
        ylim <- range(c(rsp, -rsp), na.rm = TRUE)
        if (id.n > 0) {
            ylim <- extendrange(r = ylim, f = 0.1)
            show.rsp <- order(-cook)[iid]
        }
        do.plot <- TRUE
        if (isConst.hat) {
            if (missing(caption)) 
                caption[[5]] <- "Constant Leverage:\n Residuals vs Factor Levels"
            aterms <- attributes(terms(x))
            dcl <- aterms$dataClasses[-aterms$response]
            facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
            mf <- model.frame(x)[facvars]
            if (ncol(mf) > 0) {
                effM <- mf
                for (j in seq_len(ncol(mf))) effM[, j] <- sapply(split(yh, 
                  mf[, j]), mean)[mf[, j]]
                ord <- do.call(order, effM)
                dm <- data.matrix(mf)[ord, , drop = FALSE]
                nf <- length(nlev <- unlist(unname(lapply(x$xlevels, 
                  length))))
                ff <- if (nf == 1) 
                  1
                else rev(cumprod(c(1, nlev[nf:2])))
                facval <- ((dm - 1) %*% ff)
                facval[ord] <- facval
                xx <- facval
                newplot <- xyplot(rsp ~ facval, xlim = c(-1/2, 
                  sum((nlev - 1) * ff) + 1/2), ylim = ylim, main = caption[[5]], 
                  sub = if (one.fig) 
                    sub.caption
                  else NULL, xlab = "Factor Level Combinations", 
                  ylab = ylab5, type = type, panel = function(x, 
                    y) {
                    panel.default(x, y, line.col = line.col, 
                      ...)
                    panel.abline(v = ff[1L] * (0:nlev[1L]) - 
                      1/2, col = addline.col, lty = lty[4])
                    panel.abline(h = 0, lty = lty[3], col = addline.col)
                  }, ...)
                results <- c(results, list(newplot))
            }
            else {
                message("hat values (leverages) are all = ", 
                  format(mean(r.hat)), "\n and there are no factor predictors; no plot no. 5")
                do.plot <- FALSE
            }
        }
        else {
            xx <- hii
            xx[xx >= 1] <- NA
            newplot <- xyplot(rsp ~ xx, xlim = 1.05 * c(0, max(xx, 
                na.rm = TRUE)), ylim = ylim, xlab = "Leverage", 
                ylab = ylab5, type = type, model = x, scales = list(alternating = 1, 
                  tck = c(1, 0)), main = if (one.fig) 
                  caption[[5]]
                else NULL, sub = if (one.fig) 
                  sub.caption
                else NULL, panel = function(x, y, model = x, 
                  ...) {
                  panel.abline(h = 0, v = 0, lty = lty[3], col = trellis.par.get("add.line")$col)
                  if (length(cook.levels)) {
                    p <- length(coef(model))
                    xscale <- current.viewport()$xscale
                    yscale <- current.viewport()$yscale
                    hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), 
                      max(xscale), length.out = 101)
                    for (crit in cook.levels) {
                      cl.h <- sqrt(crit * p * (1 - hh)/hh)
                      llines(hh, cl.h, lty = lty[2], col = addline.col)
                      llines(hh, -cl.h, lty = lty[2], col = addline.col)
                    }
                    grid.text("Cook's dist.", x = unit(3, "char"), 
                      y = unit(0.5, "lines"), just = c("left", 
                        "bottom"), gp = gpar(cex = 0.8, lyt = 2))
                    grid.lines(x = unit(c(0.5, 2.5), "char"), 
                      y = unit(0.7, "lines"), gp = gpar(col = addline.col, 
                        lty = lty[2], lwd = 2))
                    xmax <- min(0.99, max(xscale))
                    ymult <- sqrt(p * (1 - xmax)/xmax)
                    aty <- c(-sqrt(rev(cook.levels)) * ymult, 
                      sqrt(cook.levels) * ymult)
                    pushViewport(viewport(clip = "off", name = "unclipped5", 
                      xscale = current.viewport()$xscale, yscale = current.viewport()$yscale))
                    waty <- which(aty <= max(yscale) & aty >= 
                      min(yscale))
                    if (length(waty) > 0) {
                      grid.yaxis(at = aty[waty], label = paste(c(rev(cook.levels), 
                        cook.levels))[waty], main = FALSE, gp = gpar(cex = 0.7, 
                        tck = 0, col = addline.col))
                    }
                    upViewport()
                  }
                  if (TRUE) {
                    panel.default(x, y, type = c("p", "smooth"), 
                      line.col = line.col, ...)
                  }
                  else {
                    panel.xyplot(x, y, type = c("p"), line.col = line.col, 
                      ...)
                  }
                  if (do.plot) {
                    if (id.n > 0) {
                      grid.identify.points(xx, y, show.rsp, cex = cex.id)
                    }
                  }
                }, ...)
            results <- c(results, list(newplot))
        }
    }
    if (show[6L]) {
        g <- dropInf(hii/(1 - hii), hii)
        g <- hii/(1 - hii)
        ymx <- max(cook, na.rm = TRUE) * 1.025
        athat <- pretty(hii)
        p <- length(coef(x))
        newplot <- xyplot(cook ~ g, scales = list(alternating = 1, 
            tck = c(1, 0.3)), xlim = 1.05 * c(0, max(g, na.rm = TRUE)), 
            ylim = c(0, ymx), main = if (one.fig) 
                caption[[6]]
            else NULL, sub = if (one.fig) 
                sub.caption
            else NULL, ylab = "Cook's distance", xlab = expression("Leverage  " * 
                h[ii]), panel = function(x, y, ...) {
                pushViewport(viewport(clip = "off", name = "unclipped6", 
                  xscale = current.viewport()$xscale, yscale = current.viewport()$yscale))
                bval <- pretty(sqrt(p * cook/g), 5)
                xmax <- max(current.viewport()$xscale)
                ymax <- max(current.viewport()$yscale)
                for (i in 1L:length(bval)) {
                  bi2 <- bval[i]^2
                  if (ymax > bi2 * xmax) {
                    xi <- xmax
                    yi <- bi2 * xi
                    panel.abline(0, bi2, lty = lty[3], line.col = addline.col)
                    grid.text(paste(bval[i]), unit(1, "npc") + 
                      unit(0.3, "char"), yi, default.units = "native", 
                      just = c("left", "center"), gp = gpar(cex = 0.8, 
                        col = addline.col))
                  }
                  else {
                    yi <- ymax
                    xi <- yi/bi2
                    panel.abline(0, bi2, lty = lty[3], line.col = addline.col)
                    grid.text(paste(bval[i]), xi, unit(1, "npc") + 
                      unit(0.25, "char"), default.units = "native", 
                      just = c("center", "bottom"), gp = gpar(cex = 0.8, 
                        col = addline.col))
                  }
                }
                panel.default(x, y, line.col = line.col, ...)
                if (id.n > 0) {
                  show.r <- order(-cook)[iid]
                  grid.identify.points(g, cook, show.r, cex = cex.id)
                }
                upViewport()
            }, ...)
        results <- c(results, list(newplot))
    }
    if (print.plots) {
        for (p in results) {
            print(p)
        }
    }
    trellis.par.set(theme = old.theme)
    if (length(results) == 1) {
        invisible(results[[1]])
    }
    else {
        invisible(results)
    }
}
