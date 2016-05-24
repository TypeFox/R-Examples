## MM: more following of plot.lm() : ~/R/D/r-devel/R/src/library/stats/R/plot.lm.R

plot.lmrob <-
function (x, which = 1:5,
          caption = c("Standardized residuals vs. Robust Distances",
          "Normal Q-Q vs. Residuals", "Response vs. Fitted Values",
          "Residuals vs. Fitted Values" ,
          "Sqrt of abs(Residuals) vs. Fitted Values"),
	  panel = if(add.smooth) panel.smooth else points,
          sub.caption = deparse(x$call), main = "",
          compute.MD = TRUE, # maybe  (n < 1000 && p < 20)
          ask = prod(par("mfcol")) < length(which) && dev.interactive(),
	  id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
          label.pos = c(4,2), qqline = TRUE, add.smooth = getOption("add.smooth"),
          ..., p = 0.025)
{
    if (!inherits(x, "lmrob"))
        stop("Use only with 'lmrob' objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 5))
        stop("'which' must be in 1:5")
    show <- rep(FALSE, 5)
    show[which] <- TRUE
    r <- residuals(x)
    n <- length(r)
    sr <- r/x$scale
    yh <- fitted(x)
    if (is.null(id.n))
	id.n <- 0
    else {
	id.n <- as.integer(id.n)
	if(id.n < 0L || id.n > n)
	    stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    if(id.n > 0L) { ## label the largest residuals
	if(is.null(labels.id))
	    labels.id <- paste(1L:n)
	iid <- 1L:id.n
	show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
	## if(any(show[2L:3L]))
	##     show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
	text.id <- function(x, y, ind, adj.x = TRUE) {
	    labpos <-
		if(adj.x) label.pos[1+as.numeric(x > mean(range(x)))] else 3
	    text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
		 pos = labpos, offset = 0.25)
	}
    }

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
	if(is.null(x[['MD']]) && compute.MD) {
	    message("recomputing robust Mahalanobis distances")
	    x$MD <- ## need to recompute
		robMD(x = if(!is.null(x[['x']])) x$x else
		      if(!is.null(x[['model']])) model.matrix(x, x$model)
		      else stop("need 'model' or 'x' component for robust Mahalanobis distances"),
		      intercept = attr(x$terms,"intercept"),
                      wqr = x$qr)
	    ## try to "cache" them with the object
	    .ge <- .GlobalEnv
	    if(identical(parent.frame(), .ge) &&
	       exists((cnx <- as.character(match.call()[["x"]])), .ge)) {
		assign(cnx, x, envir = .ge)
		message("saving the robust distances 'MD' as part of ", sQuote(cnx))
	    }
	}
        if(!is.null(xD <- x[['MD']])) {
            if (p < 0 || p > 1)
                stop ("Tolerance range must be between 0% to 100%")
            else chi <- sqrt( qchisq(p = 1-p, df = x$rank) )
            ylim <- range(sr, na.rm=TRUE)
	    if(id.n > 0) ylim <- extendrange(r = ylim, f = 0.08)
	    plot(xD, xlab = "Robust Distances",
		 sr, ylab = "Robust Standardized residuals", ylim=ylim,
		 main = main, type = "n", ...)
            panel(xD, sr, ...)
            mtext(caption[1], 3, 0.25)
            if (one.fig)
                title(sub = sub.caption, ...)
	    if(id.n > 0) {
		y.id <- sr[show.r]
		y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
		text.id(xD[show.r], y.id, show.r)
	    }
            abline(h = c(2.5,-2.5), lty = 3)
            abline(v = chi, lty = 3)
        }
    }
    if (show[2L]) { ## Normal
        qq <- qqnorm(r, ylab = "Residuals", main = main,...)
        if(qqline) qqline(r, lty = 3, col = "gray50")
        mtext(caption[2], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
	if(id.n > 0)
	    text.id(qq$x[show.r], qq$y[show.r], show.r)
    }
    if (show[3]) {
        y <- if(!is.null(x[['model']])) model.response(x$model) else yh + r
        m1 <- min(yh,y)
        m2 <- max(yh,y)
        plot(yh, y, xlab = "Fitted Values", ylab = "Response",
             xlim = c(m1,m2), ylim = c(m1,m2), main = main, type = "n", ...)
        panel(yh, y, ...)
        mtext(caption[3], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
	if(id.n > 0)
	    text.id(yh[show.r], y[show.r], show.r)
        abline(a = 0,b = 1)
    }
    if (show[4]) {
        plot(yh, r, xlab = "Fitted Values", ylab = "Residuals",
             main = main, type = "n", ...)
        panel(yh, r, ...)
        mtext(caption[4], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
	if(id.n > 0) {
	    y.id <- r[show.r]
	    y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	    text.id(yh[show.r], y.id, show.r)
	}
        abline(h = c(2.5*x$scale,0,-2.5*x$scale), lty = 3)
    }
    if (show[5]) {
        sqrtabsr <- sqrt(abs(r))
        plot(yh, sqrtabsr, xlab = "Fitted Values", ylab = "Sqrt of abs(Residuals)",
             main = main, type = "n", ...)
        panel(yh, sqrtabsr, ...)
        mtext(caption[5], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
	if(id.n > 0)
	    text.id(yh[show.r], sqrtabsr[show.r], show.r)
    }
    invisible()
}

