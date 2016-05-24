n.plot <-
    function(x, y=NULL, nam = NULL, abbr = n >= 20 || max(nchar(nam))>=8,
             xlab = NULL, ylab = NULL, log = "",
             cex = par("cex"), col = par("col"), ...)
{
    ## Purpose: "Name Plot"; Names (or numbers) instead of points in plot(..)
    ## --> help(n.plot) !
    if(inherits(x,"formula")) # is(x, "formula")
        stop("formula not yet supported")
    ## this is like plot.default():
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel, log)
    xlab <- if (is.null(xlab)) xy$xlab else xlab
    ylab <- if (is.null(ylab)) xy$ylab else ylab
    plot(xy, type = 'n', xlab = xlab, ylab = ylab, log = log, ...)
    n <- length(x)
    if(is.null(nam)) {  nam <- rownames(x)
     if (is.null(nam)) { nam <- names(x)
      if (is.null(nam)) { nam <- names(y)
       if (is.null(nam)) { nam <- paste(1:n) #- Use 1,2,.. if no names
    }}}}
    if(abbr) nam <- abbreviate(nam, minlength=1)
    text(xy, labels=nam, cex=cex, col=col)
    invisible(nam)
}

TA.plot <-
  function(lm.res, fit = fitted(lm.res),
           res = residuals(lm.res, type = "pearson"),
           labels = NULL, main = mk.main(), xlab = "Fitted values",
           draw.smooth = n >= 10, show.call = TRUE, show.2sigma = TRUE,
           lo.iter = NULL, lo.cex = NULL,
           par0line  = list(lty = 2, col = "gray"),
           parSmooth = list(lwd = 1.5, lty = 4, col = 2),
           parSigma  = list(lwd = 1.2, lty = 3, col = 4),
           verbose = FALSE, ...)
{
  ## Purpose: Produce a Tukey-Anscombe plot of a linear model fit
  ##	      Note that residuals and fitted are UN-correlated (IFF intercept..)
  ## -------------------------------------------------------------------------
  ## Arguments: lm.res = RESult of lm(..)
  ##    res : (weighted) residuals by default,
  ##	labels = 'symbols' for point, default(NULL): extract names or use seq.nr
  ##             use '*' to get simple '*' symbols.
  ##
  ## --- see on-line help by  "?TA.plot" !!
  ## -------------------------------------------------------------------------
  ## Uses : n.plot(.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: Dec 92 / Nov.93;  for R: 1999/2000
  if(missing(main)) {
    call <- lm.res $ call
    if(is.call(call[["formula"]]) && any(c("lm", "aov") == call[[1]]))
      call <- call[["formula"]]
    else {  #-- only formula part; no extra  'data ='
        if (length(call) >= 3 && !is.na(m.f <- match("formula", names(call)))) {
            call <- call[c(1, m.f)]
            names(call)[2] <- ""
        }
    }
    mk.main <- function() {
      cal <- call ## if(is.R()) call else get("call", frame = sys.parent())
      if(is.null(cal))
        "Tukey-Anscombe plot of ???"
      else {
	  nc <- nchar(ccal <- deparse(cal, width.cutoff = 200)[1])
	  if(verbose)
	      cat("|cal|=", length(cal), "; nchar(ccal) =", nc,": '", ccal, "'\n", sep="")
	  if(nc > 36)
	      warning("TA.plot: 'main' title is long; consider using cex.main = 0.8",
		      call. = FALSE)
	  ##-- now should even go further:
	  ##--  E.g. if  nc > 50,  use  cex = .8 in the call to n.plot below
	  paste(if(nc < 13) "Tukey-Anscombe plot of :  "
		else if(nc < 24) "T.A. plot of: " else "TA-pl:", ccal)
	}
    }
  }
  if("ylab" %in% names(list(...))) {
      n.plot(fit, res, nam = labels, xlab = xlab, main = main, ...)
  } else {
      yl <- "Residuals"
      if(!is.null(lm.res$weights) &&
         any(abs(lm.res$resid- res) > 1e-6*mad(res)))
          yl <- paste("WEIGHTED", yl)
      n.plot(fit, res, nam = labels, xlab = xlab, ylab = yl, main = main, ...)
  }
  if(show.call)
    mtext(deparse(match.call()), side = 3, line = 0.5, cex = 0.6, adj=1)
  do.call("abline", c(list(h= 0), par0line))
  p.mgp <- par("mgp")[1:2] #-- line numbers of margin text: xlab & label
  if(missing(lo.cex))
    lo.cex <- max(.2, min(0.8*par("mex"), .9*-diff(p.mgp))/par("mfg")[4])
  m.line <- if(par("mfg")[4]==1) .2+ p.mgp[1] else
                              max(p.mgp[1] - .2*lo.cex, sum(p.mgp)/2)
  if(show.2sigma) {
    s2 <- c(-2,2) * mad(res, center=0)
    rr <- range(res)
    if(s2[1] < rr[1] || s2[2] > rr[2])
      mtext(paste("2 sigma = ", format(s2[2])),
	    side= 1, line= m.line, adj = 0, cex= lo.cex)
    ##abline(h= s2, lwd=1.8, lty=3, col=4)
    do.call("abline", c(list(h= s2), parSigma))
  }
  n <- length(res)
  if(draw.smooth) {
    if(!is.list(parSmooth)) stop("`parSmooth' must be a list")
    ##-- lo.iter: idea of Werner Stahel:  no robustness for 'glm'  residuals
    if (is.null(lo.iter))
      lo.iter <- if(inherits(lm.res, "glm")&& lm.res$family[1]!="Gaussian")
	0  else  3
    f <- max(0.2, 1.25 * n^-.2) #'-- Martin's very empirical formula...
    rlow <- lowess(fit, res, f = f, iter = lo.iter)
    do.call("lines",c(rlow, parSmooth))

    mtext(paste("-.-.-.- : lowess smooth (f =", format(round(f,2)),
		if(lo.iter!=3) paste(", it=", lo.iter), ")"),
	  side = 1, line = m.line, cex = lo.cex, adj = 1)
  }
  ##-	 "Correlation:", formatC(cor(fit,res), dig=3),
  ## mtext(paste(" -- Rank corr.:", formatC(cor(rank(fit),rank(res)), dig=3)) )
  invisible()
}
