	plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=NULL, yfun=NULL, subset=NULL,
	                       abc=NULL, add=FALSE, nlme=FALSE, ...)
{
#	plot curves from sitar model
#	opt:
#		d = fitted distance curve (labels[1] = x, labels[2] = y)
#		v = fitted velocity curve (labels[3] = y)
#		e = fitted fixed effects distance curve (labels[1] = x, labels[2] = y)
#		D = fitted distance curves by subject
#		V = fitted velocity curves by subject
#		u = unadjusted y vs t curves by subject
#		a = adjusted y vs adjusted t curves by subject
#
#		multiple options all plot on same graph
#
#		labels are for opt dv - particularly v
#		or use xlab and ylab and y2par
#
#		apv TRUE draws vertical line at age of peak velocity
#		and returns apv/pv (respecting xfun/yfun settings)
#
#		xfun/yfun are functions to apply to x/y before plotting
#
#		subset is subset of values
#
#		abc is a vector of named sitar parameters for opt dv e.g.
#			abc=c(a=1, b=0.1, c=-0.1)
#		or a single id level whose abc values are to be used
#
#		add TRUE overwrites previous graph (or use lines)
#
#		nlme TRUE plots model as nlme object

	if (nlme) {
	  do.call('plot.lme', as.list(match.call()[-1]))
	}
	else {
    options <- c('d', 'e', 'u', 'a', 'D', 'v', 'V')
    optaxis <- c(1, 1, 1, 1, 1, 2, 2) # default y axis
    optmult <- c(F, F, T, T, T, F, T) # multiple curves?
    axismin <- 3; axismax <- 0
    for (i in 1:nchar(opt)) {
      no <- match(substr(opt, i, i), options, NA)
      if (is.na(no)) next
      if (optaxis[no] > axismax) axismax <- optaxis[no]
      if (optaxis[no] < axismin) axismin <- optaxis[no]
    }
    if (axismin == 2) {
      optaxis <- optaxis - 1
      axismax <- axismin <- 1
    }
    model <- x
		data <- getData(model)
		mcall <- model$call.sitar
		x <- getCovariate(model)
		y <- getResponse(model)
		id <- getGroups(model)
		nf <- length(fitted(model))
		if (nf != length(y)) stop(paste0('model (length=', nf, ') incompatible with data (rows=', length(y), ')'))

#	extract list(...)
		ccall <- match.call()[-1]
#	subset to plot model
		subset <- eval(ccall$subset, data, parent.frame())
		if (is.null(subset)) subset <- rep(TRUE, nf)
		dots <- match.call(expand.dots=FALSE)$...
		if (length(dots) > 0) ARG <- lapply(as.list(dots), eval, data, parent.frame())
			else ARG <- NULL
#	if xlab not specified replace with label or x name (depending on xfun)
		if (!"xlab" %in% names(ARG)) {
		  if(!missing(labels)) xl <- labels[1] else {
		    if(!is.null(xfun)) xl <- paste0('(', deparse(substitute(xfun)), ')(', deparse(mcall$x), ")") else
		      xl <- ifun(mcall$x)$varname
		  }
			ARG <- c(ARG, list(xlab=xl))
		}
		else xl <- ARG$xlab
#	if ylab not specified replace with label or else y name (depending on yfun)
		if (!"ylab" %in% names(ARG)) {
		  if(!missing(labels)) yl <- labels[2] else {
		    if(!is.null(yfun)) yl <- paste0('(', deparse(substitute(yfun)), ')(', deparse(mcall$y), ")") else
		      yl <- ifun(mcall$y)$varname
		  }
			ARG <- c(ARG, list(ylab=yl))
		}
		else yl <- ARG$ylab
# if labels not specified create it
		if (missing(labels)) labels <- c(xl, yl, paste(yl, 'velocity'))
		# if (missing(labels)) labels <- c(xl, yl, ifelse(typeof(yl) == 'expression',
		#   expression(paste(as.character(yl), '~~velocity')), paste(yl, 'velocity')))

#	create output list
		xy <- list()

# derive xfun and yfun
		if (is.null(xfun)) xfun <- ifun(mcall$x)$fn
		if (is.null(yfun)) yfun <- ifun(mcall$y)$fn

#	plot y vs t by subject
		if (grepl("u", opt)) {
		  xt <- x
		  yt <- y
		  do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
		  add <- TRUE
		}

		xseq <- function(x, n=101) {
		  # n is the number of points across the x range
		  rx <- range(x, na.rm=TRUE)
		  seq(rx[1], rx[2], length.out=n)
		}

		stackage <- function(x, id, n=101) {
		  # generate x and id values across the x range to plot spline curves
		  npt <- n / diff(range(x))
		  xid <- by(data.frame(x=x, id=id), id, function(z) {
		    nt <- floor(npt * diff(range(z$x))) + 1
		    data.frame(x=seq(min(z$x), to=max(z$x), length.out=nt), id=rep.int(z$id[[1]], nt))
		  })
		  df <- xid[[1]][FALSE, ]
		  for (dft in xid) df <- rbind(df, dft)
		  df
		}

#	plot fitted curves by subject
		if (grepl("D", opt)) {
		  newdata=stackage(x[subset], id[subset])
		  newdata <- cbind(newdata, y=predict(model, newdata=newdata, xfun=xfun, yfun=yfun))
		  do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
		                          data=newdata, add=add), ARG))
		  add <- TRUE
		}

#	plot fitted velocity curves by subject
		if (grepl("V", opt)) {
		  newdata=stackage(x[subset], id[subset])
		  newdata <- cbind(newdata, y=predict(model, newdata=newdata, deriv=1, xfun=xfun, yfun=yfun))
		  ARG$ylab <- labels[3]
		  do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
                              data=newdata, add=add), ARG))
      add <- TRUE
		}

#	plot fitted distance and velocity curves
		if (grepl("d", opt) || grepl("v", opt) || apv) {
			xt <- xseq(x[subset])
  		newdata <- data.frame(x=xt)
# if subset, flag for predict
      if (!identical(subset, rep(TRUE, nf))) attr(newdata, 'subset') <- subset

#	adjust for abc
			if (!is.null(abc)) {
#	abc is id level
				if (length(abc) == 1) {
					idabc <- rownames(ranef(model)) %in% abc
					if (sum(idabc) == 0) stop(paste('id', abc, 'not found'))
					abc <- ranef(model)[idabc, ]
				}
#	abc is named vector
			  else if (length(abc) > 3 || is.null(names(abc)))
			    stop('abc should be either single id level or up to three named random effect values')
			}
  		else abc <- ranef(model)

			yt <- yfun(predict(object=model, newdata=newdata, level=0, abc=abc))
			vt <- predict(object=model, newdata=newdata, level=0, deriv=1, abc=abc, xfun=xfun, yfun=yfun)
#	derive cubic smoothing spline curve
			xt <- xfun(xt)
			xy$ss <- ss <- makess(xt, yt)

			if (grepl("d", opt) && grepl("v", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, y2=vt, labels=labels, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("d", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("v", opt)) {
				ARG$ylab <- labels[3]
				xy <- do.call("y2plot", c(list(x=xt, y1=vt, labels=labels[c(1,3)], add=add, xy=xy), ARG))
				add <- TRUE
			}
		}

#	plot fixed effects distance curve
		if (grepl("e", opt)) {
  		xt <- xseq(x[subset])
      yt <- predict(model$ns, newdata=data.frame(x=xt - model$xoffset))
			ox <- order(xt)
			xy <- do.call("y2plot", c(list(x=xfun(xt[ox]), y1=yfun(yt[ox]), add=add, xy=xy), ARG))
			add <- TRUE
		}

#	plot adjusted y vs adjusted t by subject
		if (grepl("a", opt)) {
			yt <- xyadj(x, y, id, model)
			xt <- yt$x
			yt <- yt$y
    	do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
			add <- TRUE
		}

#	plot vertical line at age of peak velocity
		if (apv) {
			xy$apv <- ss$apv
			if (!is.na(opt)) print(signif(xy$apv, 4))
			if (add) {
				if (is.null(ARG$y2par$lty)) ARG$y2par$lty <- 3
				do.call('abline', c(list(v=xy$apv[1]), ARG$y2par))
			}
		}
		invisible(xy)
	}
}
