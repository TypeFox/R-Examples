loessDemo <-
    function(x, y, span = 1/2, degree = 1, family = c("gaussian", "symmetric"),
	     nearest = FALSE, nout = 501,
	     xlim = numeric(0), ylim = numeric(0), strictlim=TRUE, verbose = TRUE,
	     inch.sym = 0.25, pch = 4,
	     shade = TRUE, w.symbols = TRUE,
	     sym.col = "blue", w.col = "light blue", line.col = "steelblue")
{
    ## function to demonstrate the locally weighted regression function loess
    ## written and posted to S-news, Thu, 27 Sep 2001 07:48
    ###		Dr. Greg Snow
    ## 		Brigham Young University, Department of Statistics
    ## 		gls@byu.edu
    ## Modified by Henrik Aa. Nielsen, IMM, DTU (han@imm.dtu.dk)
    ## spiffed up (R only), by M.Mächler, SfS ETH Zurich

    family <- match.arg(family)
    ## drop NA's and sort:
    miss.xy <- is.na(x) | is.na(y)
    x <- x[!miss.xy]
    y <- y[!miss.xy]
    ix <- order(x)
    x <- x[ix]
    y <- y[ix]
    degree <- as.integer(degree)
    if(length(degree) != 1 || is.na(degree) || degree < 0 || 2 < degree)
	stop("'degree' must be in {0,1,2}")

    fit.D <- loess(y ~ x, degree = degree, span = span, family = family,
		   control = loess.control(surface = "direct"))

    fit.I <- loess(y ~ x, degree = degree, span = span, family = family)

    xx <- seq(min(x), max(x), len = nout)
    est <- list(x = xx,
		y = predict(fit.I, newdata = data.frame(x = xx)))

    xl <- if(strictlim && is.numeric(xlim) && length(xlim) == 2)
	      xlim
	  else {
	      xl <- range(x, est$x, xlim)
	      xl <- xl + c(-1, 1) * 0.03 * diff(xl)
	  }
    yl <- if(strictlim && is.numeric(ylim) && length(ylim) == 2) {
	      dy <- 0.05 * diff(ylim)
	      ylim
	  }
	  else {
	      yl <- range(y, est$y, ylim, fitted(fit.D))
	      dy <- 0.05 * diff(yl)
	      yl + c(-1, 1) * dy
	  }
    ## room below for weights
    dy <- 4*dy
    yl[1] <- yl[1] - dy
    stit <- paste("span = ", span,";  degree = ", degree)
    if(family != "gaussian")
	stit <- paste(stit,".  family = \"", family,'"',sep="")

    fitPlot <- function(x, y, w, est, fit.D, xl, yl)
    {
	pU <- par("usr")
	plot(x, y, pch = pch, xlim = xl, ylim = yl, sub = stit)
	if(!is.null(w)) {
	    w <- w/max(w) # in [0,1]
	    wP <- w > 1e-5
	    nw <- length(xw <- x[wP])
	    if(w.symbols)
		symbols(xw, y[wP], circles = sqrt(w[wP]),
			inches = inch.sym, add = TRUE, fg = sym.col)
	    # scale [0,1] to yl[1] + [0, dy] :
	    y0 <- pU[3]
	    wy <- y0 + (dy+yl[1]-y0) * w[wP]
	    polygon(c(xw[1], xw, xw[nw]), c(y0, wy, y0), col = w.col)
	    segments(xw, rep(y0,nw), xw, wy, col=sym.col)
	}
	lines(x, fitted(fit.D), col = 2, lwd = 2)
	mtext("Exact estimate with linear interpolation between x-values ('surface = \"direct\")",
	      col = 2, adj = 0.5, line = 0.5)
	lines(est, col = 3, lwd = 2)
	mtext("Estimate obtained using the default interpolation scheme",
	      col = 3, adj = 0.5, line = 2)
	pU
    }

    fitPlot(x, y, w=NULL, est, fit.D, xl, yl)

     repeat {
	if(verbose)
	    cat("click left for x0  to predict -- click right to stop ")
	x0 <- locator(1)$x
	if(verbose) cat("\n")
	if(length(x0) < 1)## right clicking leaves loop
	    break
	if(nearest)
	    x0 <- unique(x[abs(x - x0) == min(abs(x - x0))])
	if(verbose)
	    cat("x0 =", x0, "\n")
	Dx <- abs(x - x0)
	d <-
	    if(span < 1)
		sort(Dx)[as.integer(span * length(x))]
	    else max(Dx) * sqrt(span)
	w <- rep(0, length(x))
	s <- Dx <= d
	w[s] <- (1 - (Dx[s]/d)^3)^3 # tricube weights
	pU <- fitPlot(x, y, w, est, fit.D, xl, yl)
	##    =======                  ==

	if(degree > 0L) { ## is '1' or '2
	    if(degree == 1)
		abline(lm(y ~ x, weights = w), col = line.col)
	    else ## (degree == 2) # predict(lm( ~ poly()) fails!
		lines(xx, predict(lm(y ~ x + I(x^2), weights = w),
				  data.frame(x=xx)),
		      col = line.col, err = -1)
	} else { ## degree == 0
	    ##lines(x, fitted(lm(y ~ 1, weights = w)), col = line.col, err = -1)
	    abline(a = sum(w*y)/sum(w), b = 0, col = line.col)
	}
	abline(v = x0, col = line.col, lty = 3, lwd = 0.2)
	axis(1, at= x0, labels = formatC(x0, digits=3), col.axis = line.col)
	if((x1 <- x0 - d) > xl[1]) {
	    abline(v = x1, col = line.col, lty = 2)
	    if(shade) polygon(c(pU[1],x1,x1,pU[1]), pU[c(3,3, 4,4)], density = 5)
	}
	if((x1 <- x0 + d) < xl[2]) {
	    abline(v = x1, col = line.col, lty = 2)
	    if(shade)
		polygon(c(x1, pU[c(2,2)],x1), pU[c(3,3, 4,4)], density = 5,
			angle = -45)
	}
    }
}
