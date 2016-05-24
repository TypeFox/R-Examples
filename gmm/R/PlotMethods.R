#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

plot.gmm <- function (x, which = c(1L:3),
	    main = list("Residuals vs Fitted values", "Normal Q-Q",
	    "Response variable and fitted values"),
	    panel = if(add.smooth) panel.smooth else points,
	    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
	    add.smooth = getOption("add.smooth"))
{
    if (!inherits(x, "gmm"))
	stop("use only with \"gmm\" objects")
    if (!inherits(x, "gel"))
	{
    	if(!is.numeric(which) || any(which < 1) || any(which > 3))
		stop("'which' must be in 1L:3")
	show <- rep(FALSE, 3)
	}
    else
	{
	if(!is.numeric(which) || any(which < 1) || any(which > 4))
		stop("'which' must be in 1L:4")
	show <- rep(FALSE, 4)
	}
    
    show[which] <- TRUE
    r <- residuals(x)
    if(ncol(r)>1)
	stop("plot.gmm is not yet implemented for system of equations")

    yh <- fitted(x) 
    n <- length(r)
    
    if (ask) {
	oask <- devAskNewPage(TRUE)
	on.exit(devAskNewPage(oask))
    }

    ##---------- Do the individual plots : ----------

    if (show[1L]) {
	ylim <- range(r, na.rm=TRUE)
	ylim <- extendrange(r= ylim, f = 0.08)
	plot(yh, r, xlab = "Fitted", ylab = "Residuals", main = main[1L],
	     ylim = ylim, type = "n", ...)
	panel(yh, r)
	abline(h = 0, lty = 3, col = "gray")
    }
    if (show[2L]) { ## Normal
	rs <- (r-mean(r))/sd(r)
	ylim <- range(rs, na.rm=TRUE)
	ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
	qq <- qqnorm(rs, main = main[2L], ylab = "stand. residuals", ylim = ylim, ...)
	qqline(rs, lty = 3, col = "gray50")
    }
    if (show[3L]) {
	y <- as.matrix(model.response(x$model, "numeric"))
	ylim <- range(yh, na.rm=TRUE)
	ylim <- extendrange(r= ylim, f = 0.08)
	plot(y, main = main[3L],
	     ylim = ylim,  ...)
	lines(yh,col=2)
    }
    if (inherits(x, "gel"))
	{
	    if (show[4L]) {
		pt <- x$pt
		plot(pt, type='l',main = main[4L],ylab="Implied Prob.", ...)
		emp_pt <- rep(1/length(pt),length(pt))
		lines(emp_pt,col=2)
		legend("topleft",c("Imp. Prob.","Empirical (1/T)"),col=1:2,lty=c(1,1))
    		}
	}


    invisible()
}

plot.gel <- function (x, which = c(1L:4),
	    main = list("Residuals vs Fitted values", "Normal Q-Q",
	    "Response variable and fitted values","Implied probabilities"),
	    panel = if(add.smooth) panel.smooth else points,
	    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
	    add.smooth = getOption("add.smooth"))
	{
	if (!inherits(x, "gel"))
		stop("use only with \"gel\" objects")	
	class(x) <- c("gmm","gel")
	plot(x,which=which,main=main,panel=panel, ask=ask, ..., add.smooth=add.smooth)
	
	invisible()
	}
