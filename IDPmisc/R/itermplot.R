## itermplot.R


itermplot <- function(model, data=NULL,envir=environment(formula(model)),
                      partial.resid=FALSE, scale=0,
                      pixs=1, zmax=NULL, ztransf = function(x){x},
                      colramp = IDPcolorRamp,
                      terms=NULL, se=FALSE,
                      xlabs=NULL, ylabs=NULL,
                      main = NULL, col.term = "black", lwd.term = 2,
                      col.se = "gray", lty.se = 2, lwd.se = 1,
                      col.smth = "darkred", lty.smth=2, lwd.smth=2,
                      span.smth=2/3,
                      ask = interactive() && nb.fig < n.tms &&
                      .Device != "postscript",
                      use.factor.levels=TRUE, smooth=NULL,
                      ...)
  ## modified Version of termplot R2.3.1
  ## Author:  Rene Locher
  ## Version: 2007-08-14
{
    which.terms<-terms
    terms <- ## need if(), since predict.coxph() has non-NULL default terms :
	if (is.null(terms))
	    predict(model, type="terms", se.fit=se)
	else
	    predict(model, type="terms", se.fit=se, terms=terms)
    n.tms <- ncol(tms <- as.matrix(if(se) terms$fit else terms))
    mf <- model.frame(model)
    if (is.null(data))
        data<-eval(model$call$data,envir)
    if (is.null(data))
        data<-mf
    if (NROW(tms)<NROW(data)){
        use.rows<-match(rownames(tms),rownames(data))
      } else use.rows<-NULL
    nmt <- colnames(tms)
    cn <- parse(text=nmt)
    ## Defaults:
    if (!is.null(smooth))
      smooth<-match.fun(smooth)
    if (is.null(ylabs))
	ylabs <- paste("Partial for",nmt)
    if (is.null(main))
        main <- ""
    else if(is.logical(main))
        main <- if(main) deparse(model$call, 500) else ""
    else if(!is.character(main))
        stop("'main' must be TRUE, FALSE, NULL or character (vector).")
    main <- rep(main, length.out = n.tms) # recycling
    pf <- envir
    carrier <- function(term) { # used for non-factor ones
	if (length(term) > 1)
	    carrier(term[[2]])
	else
	    eval(term, data, enclos = pf)
    }
    carrier.name<-function(term){
      	if (length(term) > 1)
	    carrier.name(term[[2]])
	else
	    as.character(term)
    }
    if (is.null(xlabs))
        xlabs<-unlist(lapply(cn,carrier.name))

    if (partial.resid || !is.null(smooth)){
	pres <- residuals(model, "partial")
        if (!is.null(which.terms)) pres<-pres[,which.terms,drop=FALSE]
      }
    is.fac <- sapply(nmt, function(i) is.factor(mf[,i]))

    se.lines <- function(x, iy, i, ff = 2) {
        tt <- ff * terms$se.fit[iy,i]
        lines(x, tms[iy,i] + tt, lty=lty.se, lwd=lwd.se, col=col.se)
        lines(x, tms[iy,i] - tt, lty=lty.se, lwd=lwd.se, col=col.se)
    }

    nb.fig <- prod(par("mfcol"))
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }

    cntsmax <- 0 ## maximum number of points per pixel

    ##---------- Do the individual plots : ----------

    for (i in 1:n.tms) {
	ylims <- range(tms[,i], na.rm=TRUE)
	if (se)
	    ylims <- range(ylims,
			   tms[,i] + 1.05*2*terms$se.fit[,i],
			   tms[,i] - 1.05*2*terms$se.fit[,i], na.rm=TRUE)
	if (partial.resid)
	    ylims <- range(ylims, pres[,i], na.rm=TRUE)
        if (scale>diff(ylims))
            ylims <- ylims+c(-0.5,0.5)*(scale-diff(ylims))
	if (is.fac[i]) {
	    ff <- mf[,nmt[i]]
            if (!is.null(model$na.action))
              ff<-naresid(model$na.action,ff)
	    ll <- levels(ff)
	    xlims <- range(seq(along.with=ll)) + c(-.5, .5)
            xx <- as.numeric(ff) ##need if rug or partial

	    plot(1,0, type = "n", xlab = xlabs[i], ylab = ylabs[i],
                 xlim = xlims, ylim = ylims, main = main[i],xaxt="n", ...)
            if (use.factor.levels)
                axis(1,at=seq(along.with=ll),labels=ll,...)
            else
                axis(1)
            if (partial.resid){
              if (!is.fac[i] && !is.null(smooth)){
                cntsmax <-
                  max(cntsmax,
                      smooth(xx,pres[,i],
                             lty=lty.smth, col=col.smth, lwd=lwd.smth,
                             span=span.smth, pixs = pixs, zmax = zmax,
                             ztransf = ztransf, colramp =colramp))
              } else
             ## points(xx, pres[,i], cex = cex, pch = pch, col = col.res)
              cntsmax <-
                max(cntsmax,
                    Image(xx, pres[,i], pixs = pixs, zmax = zmax,
                          ztransf = ztransf, colramp = colramp,
                          factors = c(is.fac[i],FALSE)))
            }
	    for(j in seq(along.with=ll)) {
		ww <- which(ff==ll[j])[c(1,1)]
		jf <- j + c(-.4, .4)
		lines(jf,tms[ww,i], col=col.term, lwd=lwd.term, ...)
		if(se) se.lines(jf, iy=ww, i=i)
	    }
	}
	else { ## continuous carrier
	    xx <- carrier(cn[[i]])
            if (!is.null(use.rows)) xx<-xx[use.rows]
	    xlims <- range(xx,na.rm=TRUE)
	    oo <- order(xx)
	    plot(xx[oo], tms[oo,i], type = "n", xlab = xlabs[i],
                 ylab = ylabs[i],
		 xlim = xlims, ylim = ylims, main = main[i])
            if (partial.resid){
              if (!is.fac[i] && !is.null(smooth)){
                cntsmax <-
                  max(cntsmax,
                      smooth(xx,pres[,i],
                             lty=lty.smth, col=col.smth, lwd=lwd.smth,
                             span=span.smth, pixs = pixs, zmax = zmax,
                             ztransf = ztransf, colramp =colramp))
              } else
             ## points(xx, pres[,i], cex = cex, pch = pch, col = col.res)
              cntsmax <-
                max(cntsmax,
                    Image(xx, pres[,i], pixs = pixs, zmax = zmax,
                          ztransf = ztransf, colramp = colramp,
                          factors = c(is.fac[i],FALSE)))
            }
            lines(xx[oo], tms[oo,i], col=col.term, lwd=lwd.term, ...)
            if(se) se.lines(xx[oo], iy=oo, i=i)
          }

      }
    invisible(zmax)
  } ## itermplot
