## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' @title Creating Nice RGB Colors
##' @param n number of colors
##' @param cols hex values or names of colors
##' @return hex values of RGB colors
##' @author Marius Hofert
##' @note not color blind safe
myRGB <- function(n, cols=c("#FFB686", "#F083AB", "#BB65C4", "#0066CC")) # myHCL(4)
    colorRampPalette(cols, space="Lab")(n) # use nice Lab space

##' @title Creating Nice HCL Colors
##' @param n number of colors
##' @param cols hex values of start and end colors
##' @return hex values of HCL colors
##' @author Marius Hofert
myHCL <- function(n, cols=c("#FFB686", "#0066CC")) {
    stopifnot(length(cols)==2)
    start <- rev(as(hex2RGB(cols[1]), "polarLUV")@coords) # start color
    end <- rev(as(hex2RGB(cols[2]), "polarLUV")@coords) # end color
    end[1] <- end[1]-360 # map correctly to not run over green
    heat_hcl(n, h=c(start[1], end[1]),
	       c.=c(start[2], end[2]),
		l=c(start[3], end[3])) # all colors on the path
}

##' for "plot math" labels {from dimnames etc}:
tryExpr <- function(ch) {
    if(inherits(r <- tryCatch( parse(text=ch, srcfile=NULL),
			      error=function(e)e),
		"error"))
	ch else r
}

##' @title Compute a grob for a background with grid lines
##' @param v x-axis locations for the vertical grid lines
##' @param h y-axis locations for the horizontal grid lines
##' @param default.units default units
##' @param name name of the grob/grobTree
##' @param gp graphical parameters
##' @param vp viewport
##' @return background-with-grid-lines grob
##' @author Marius Hofert
bgGrob <- function(v, h, default.units="native", name=NULL, gp=gpar(), vp=NULL)
    grobTree(rectGrob(),
             segmentsGrob(v, unit(0, "npc"), v, unit(1, "npc"),
                          default.units=default.units),
             segmentsGrob(unit(0, "npc"), h, unit(1, "npc"), h,
                          default.units=default.units),
             name=name, gp=gp, vp=vp)

##' @title Function for constructing a \dQuote{pair} of panel function and corresponding legend grob
##' @param method panel function method (boxplot, lines, points,...); currently must be one of
##'  "boxplot", "lines".  Instead of "points", use "lines" with \code{type="b"}
##' @param type see lines() (for the panel function)
##' @param pch see legendGrob()
##' @param lwd line width for panel
##' @param cex (plot) character size (expansion) for panel
##' @param do.legend logical indicating if a legend should be added ....{FIXME}
##' @param leg.lwd line width for legend symbols
##' @param leg.cex size (expansion) for legend symbols and labels
##' @param leg.col legend colors
##' @param leg.expr legend labels as expressions
##' @param verbose option for marking NAs
##' @param ... additional parameters passed to the panel function
##' @author Marius Hofert
panelLegend <- function(method=c("boxplot", "lines"), type, pch,
                        lwd = 1, cex = 1, do.legend=TRUE,
                        leg.lwd = lwd, leg.cex = cex, leg.col, leg.expr,
                        panel.first = NULL, panel.last = NULL,
                        verbose=getOption("verbose"), ...)
{
    ## for marking NAs
    markNA <- function(x.na, col) {
        axis(1, at=x.na, labels=rep.int("N", length(x.na)),
             mgp = c(1.5, 0.6, 0),
             tick=FALSE, cex.axis=0.9, col.axis=adjustcolor(col, 0.8))
        if(verbose)
            message("NA's at x = ", paste(format(x.na, collapse=", ")))
    }

    ## panel function
    method <- match.arg(method)
    panel <- switch(method,
                    "boxplot" = {
                        function(x, y, col, ...) {
                            stopifnot(is.matrix(y), (p <- length(x)) == ncol(y))
			    wex <- 0.5 * diff(range(x))/ p
                            if (!is.null(panel.first)) panel.first(x, y, col, ...)
			    boxplot.matrix(y, at = x, boxwex = wex, border = col,
					   col = NA, axes = FALSE, add = TRUE, names = FALSE,
					   ## cex, lwd, pch: "globals" passed here:
					   cex = cex, lwd = lwd, pch = pch, ...)
                            my <- apply(y, 2, median, na.rm=TRUE)
                            lines(x, my, col=col, lwd=lwd) # lines
                            if(any(iy <- is.na(y))) { ## Mark NA's at the bottom
                                markNA(x[apply(iy, 2, any)], col)
                            }
                            if (!is.null(panel.last)) panel.last(x, y, col, ...)
                        }
                    },
                    "lines" = {
                        function(x, y, col, ...) {
                            ## lines/points
                            if (!is.null(panel.first)) panel.first(x, y, col, ...)
                            lines(x, y, type=type, pch=pch, col=col, lwd=lwd, cex=cex, ...)
                            if(any(iy <- is.na(y))) markNA(x[iy], col) # Mark NA's at the bottom
                            if (!is.null(panel.last)) panel.last(x, y, col, ...)
                        }
                    },
                    stop(gettextf("method '%s' not yet implemented", method), domain=NA))

    ## legend function
    lg <- if(do.legend) {
        legendGrob(as.expression(leg.expr), nrow=1, pch=pch,
                   gp=gpar(col=leg.col, lwd=lwd, cex=cex))
        ## => if lwd is omitted, do.lines = FALSE => *nothing* plotted!
    }
    else nullGrob()

    ## return
    list(panel=panel, legend=lg)
}
##' { end } panelLegend

##' Matrix plot with grid and gridBase, see help(mayplot) ../man/mayplot.Rd
##' { mayplot }
mayplot <- function(x, vList, row.vars = NULL, col.vars = NULL, xvar,
                    method = if(has.n.sim) "boxplot" else "lines",
                    panel.first = NULL, panel.last = NULL,
                    type = "l", pch = NULL,  ylim = "global",
                    log = "", do.legend = TRUE,
		    spc = c(0.04/max(1,n.x-1), 0.04/max(1,n.y-1)),
                    axlabspc = c(0.12, 0.08), labspc = c(0.04, 0.04),
                    n.sim.spc = 0.06,
                    auxcol = c("gray40", "gray78", "gray90", "white"),
		    pcol = c("black", "blue", "red", "orange"), grid.lwd=1.6, ax.lwd=2,
		    tx.cex = 1.2, leg.cex = 1, xlab = NULL, ylab = NA,
                    do.n.sim = has.n.sim,
		    verbose = getOption("verbose"), show.layout = verbose, ...)
{
    ## basics ##################################################################

    avar <- names( dn <- dimnames(x) ) # *a*ll variables
    if(!all(avar %in% names(vList)))
	stop("names(dimnames(.)) should all appear in names(vList)")
    d <- dim(x)
    if(!identical(avar, nd <- names(d))) {
	if(!is.null(nd)) {
	    nonB <- nd != ""
	    if(any(nd[nonB] != avar[nonB]))
		warning("names(dim(x)) do not match names(dimnames(x))")
	}
	names(d) <- avar
    }

    ## deal with n.sim
    has.n.sim <- "n.sim" %in% avar
    if(has.n.sim)	 {
	vLnsim <- vList[["n.sim"]]
	if(vLnsim[["type"]] != "N")
	    message("\"n.sim\"'s \"type\" is not \"N\".  Hopefully on purpose")
	stopifnot((n.sim <- d[["n.sim"]]) == vLnsim[["value"]])
    }
    else if(method == "boxplot")
	stop("method 'boxplot' not possible with no 'n.sim'")

    ## deal with boxplot case
    if(method == "boxplot" && n.sim <= 4)
	warning("box plots with ", n.sim, " values are rarely useful")
    else if(method != "boxplot" && has.n.sim && n.sim > 1)
	stop("has n.sim (= ", n.sim, ") > 1, but no boxplots.",
	     "\n *** error bars not yet implemented!")

    ## variables displayed in a single panel (different colors)
    ovar <- c(xvar, col.vars, row.vars, if(has.n.sim) "n.sim") # *o*uter variables
    if(!all(ovar %in% avar)) { ## error
	if(!all(i <- xvar %in% avar))
	    stop(gettextf("'%s' entry *not* in names(dimnames(x)): \"%s\"\n",
			  "xvar", xvar[which(!i)[1]]), domain = NA)
	if(!all(i <- col.vars %in% avar))
	    stop(gettextf("'%s' entry *not* in names(dimnames(x)): \"%s\"\n",
			  "col.vars", col.vars[which(!i)[1]]), domain = NA)
	if(!all(i <- row.vars %in% avar))
	    stop(gettextf("'%s' entry *not* in names(dimnames(x)): \"%s\"\n",
			  "row.vars", row.vars[which(!i)[1]]), domain = NA)
    }
    v <- setdiff(avar, ovar) # variables displayed in a single panel
    lv <- length(v) # number of different variables ('n', 'd',...) in a single panel
    if(lv == 0) {
        if(verbose) cat("no 'inner plot' variable 'v'\n")
	v <- NULL
	v.col <- pcol[1]
        nv <- 0
    } else if(lv == 1) {
        nv <- length(dn[[v]]) # number of values the (single) panel plot variable 'v' takes
        stopifnot(nv >= 1)
	v.col <- colorRampPalette(pcol, space="Lab")(nv) # colors for v
    } else stop("length(v)=", lv, " not yet implemented; too many variables (per plot panel)")
    ## TODO MM: great if this worked with length(v) == 2 {"n" & "d"}, and maybe > 2
    ##          (e.g. *overplot*: one full plot for each "n" )
    ##          => need to have accordingly many legend lines then

    ## permute so that we know the indices of v and outer variables (may be NULL)
    avar.sort <- c(row.vars, col.vars, v, if(has.n.sim) "n.sim", xvar)
    x <- aperm(x, perm = avar.sort) # permute x

    ## verbose
    if(verbose) {
	cat("xvar, col.vars, row.vars, v, n.sim : those *not* NULL are:\n")
	print(c(xvar=xvar, col.vars=col.vars, row.vars=row.vars, v=v, n.sim=if(has.n.sim) n.sim))
	cat("method = ", method, "; ", sep="")
	stopifnot(identical(dim(x), unname(d[avar.sort])))
	cat("dim(x {after aperm()}):\n"); print(d[avar.sort])
    }

    ## z = value of xvar or seq_along() of it
    xv <- vList[[xvar]]$value # *x* axis *v*alue
    z <- if(is.character(xv) || is.recursive(xv)) seq_along(xv)
    else {
	z <- tryCatch(as.numeric(xv), error=function(e) e)
	if(inherits(z, "error")) seq_along(xv) else z
    }
    zrange <- range(z) # for forcing the same x axis limits per row
    if(method=="boxplot") zrange <- extendrange(zrange)

    ## ylim
    ylim.msg <- "'ylim' must be \"global\", \"local\", or numeric of length two"
    if(is.character(ylim)) {
	ylim <- match.arg(ylim, c("global", "local"))
	do.ylim <- ylim == "local"
	if(!do.ylim) ylim <- range(x, na.rm=TRUE)
    } else if(is.numeric(ylim)) {
	do.ylim <- FALSE
	if(!(length(ylim) == 2 && is.finite(ylim)))
	    stop(ylim.msg)
    } else stop(ylim.msg)

    do.legend <- do.legend && lv
    ## number of column/row variables
    getdn <- function(nm) if(is.null(nm)) nm else dn[[nm]] # get dim names
    nx <- length(getdn(col.vars)) # number of column variables
    ny <- length(getdn(row.vars)) # number of row variables
    ## adjust for degenerate case
    n.x <- max(1, nx) # number of columns to plot
    n.y <- max(1, ny) # number of rows to plot

    ## panel (function) *and* legend (grob) ####################################
    check.panel <- function(p)
        is.function(p) && length(ff <- formals(p)) >= 3 &&
            c("y", "col") %in% names(ff[-1]) # all but first must match

    ## construct default panel function with matching legend
    if(is.null(pch)) pch <- if(type %in% c("p","o","b")) 1 # else NULL; pch from legendGrob()
    ## note: nv can be 0 in which case lexpr == list() [still works]
    lexpr <- lapply(seq_len(nv), function(k)
                    substitute(expr == val,
                               list(expr=vList[[v]][["expr"]],
                                    val=dn[[v]][k]))) # legend symbols and labels
    pl <- panelLegend(method=method, type=type, pch=pch,
                      do.legend=do.legend,
                      leg.col=v.col, leg.expr=lexpr,
                      panel.first=panel.first, panel.last=panel.last,
                      verbose=verbose, leg.cex=leg.cex, ...)
    ## panel
    panel <- pl$panel
    stopifnot(check.panel(panel))
    ## legend
    if(do.legend) {
        legend.grob <- pl$legend
        legspc <- 0.03 ### FIXME from grob! (and test with cex=0.5 or 2,etc)
    }

    stopifnot(length(spc) == 2, length(axlabspc) == 2,
	      length(labspc) == 2, length(auxcol) == 4)

    ## number of columns/rows in the layout
    ## for y axis label; (nx-1)* (plot+gaps); plot; row.vars labels :
    ## space for 1 (of n.x) plot in x direction:
    w.all <- (axlabspc[1] + labspc[1] + n.sim.spc + (n.x-1)*spc[1])
    if(w.all >= 1) stop("Choose smaller 'axlabspc[1]', 'labspc[1]', 'spc[1]', 'n.sim.spc'")
    Pxsp <- (1 - w.all) / n.x
    wid.gl <- c(axlabspc[1], # y axis label
                rep(c(Pxsp, spc[1]), n.x-1), Pxsp, # panels (+ horizontal space in between)
                labspc[1], # panel labels
                if(has.n.sim) n.sim.spc) # if n.sim, 'second y axis' label
    nx. <- length(wid.gl) # number of columns in the layout
    ## for col.vars labels; (ny-1)*(plot+gap); plot;
    ## space for 1 (of n.y) plot in y direction:
    leg.spc <- if(do.legend) legspc else 0
    h.all <- (labspc[2] + axlabspc[2] + leg.spc + (n.y-1)*spc[2])
    if(h.all >= 1) stop("Choose smaller 'xlabspc[2]', 'labspc[2]', 'spc[2]', 'legspc'")
    Pysp <- (1 - h.all) / n.y
    hei.gl <- c(labspc[2], # panel labels
                rep(c(Pysp, spc[2]), n.y-1), Pysp, # panels (+ vertical space in between)
                axlabspc[2], # x axis label
		if(do.legend) leg.spc) # legend
    ny. <- length(hei.gl) # number of rows in the layout


    ## layout via grid #########################################################

    ## plot settings, restored on exit
    opar <- par(no.readonly=TRUE); on.exit(par(opar))
    plot.new() # start (empty) new page with 'graphics'; grid.newpage() fails here
    ## TODO MM: use smart units for widths/heights rather than just 'npc'
    gl <- grid.layout(nrow=ny., ncol=nx., widths=wid.gl, heights=hei.gl,
                      ## units in npc as for pdf()
                      ## (no square plotting region otherwise)
                      default.units="npc")
    if(show.layout) grid.show.layout(gl, vp=viewport(width=1.25, height=1.25))
    pushViewport(viewport(layout=gl)) # use this layout in a viewport
    ## drawback of working with base graphics: layout has to be fully
    ##                                         specified before plotting


    ## plot data ###############################################################

    ## walk over rows
    for(i in 1:n.y)
    {
        i. <- 2*i # column index in layout (for jumping over gaps)
	if(verbose) cat(sprintf("plot row %d (%d): [columns:] ", i, i.))
        ## pick out data to plot
	xi <- # x for a fixed 'i'
	    if(ny) switch(length(d),
			  x[i  ,drop=FALSE],
			  x[i , ,drop=FALSE],
			  x[i ,, ,drop=FALSE],
			  x[i ,,, ,drop=FALSE],
			  x[i ,,,, ,drop=FALSE],
			  stop(gettextf("'x' has wrong rank = length(dim(x)) = %d",
                                        length(d))))
	    else x
	if(do.ylim) ylim <- range(xi, na.rm=TRUE) # forcing same y axis limits per row

        ## walk over columns
	for(j in 1:n.x)
        {
	    j. <- 2*j # row index in layout (for jumping over gaps)
	    if(verbose) cat(sprintf("%d (%d) ", j, j.))

            ## push viewport
            pushViewport(viewport(layout.pos.row=i., layout.pos.col=j.))

            ## plot panels #####################################################

            ## start a base graphics plot
	    par(plt=gridPLT(), new=TRUE) # Paul: always do this before each new base graphics plot
	    ## set up coordinate system
	    plot.window(zrange, ylim, log=log)

            ## setup for background grid
            ## former base-graphics grid:
            ##   grid.rect(gp=gpar(col=NA, fill=auxcol[3]))
	    ##   grid(col=auxcol[4], lty="solid", lwd=grid.lwd, equilogs=FALSE)
            pxlog <- par("xlog")
            pylog <- par("ylog")
            v <- axTicks(1, axp=par("xaxp"), log=pxlog) # x values of vertical lines
            h <- axTicks(2, axp=par("yaxp"), log=pylog) # y values of horizontal lines
            if(pxlog) v <- log10(v)
            if(pylog) h <- log10(h)
            ## => v, h have to be determined *here* [not in bgGrob()]
            usr <- par("usr")

            ## background
            grid.draw(bgGrob(v=v, h=h, vp=viewport(xscale=usr[1:2], yscale=usr[3:4]),
                             gp=gpar(col=auxcol[4], fill=auxcol[3])))

            ## panel
	    nxy <- nx && ny # TRUE <=> there is at least one row and one column
	    if(lv) {
                for(k in 1:nv) {# walk over all variables plotted in a single panel
                    ## y values
                    y. <- switch(length(d)-2, # length(ovar)-1 - has.n.sim,
                                 if(nxy) xi[1,j,k]   else if(ny) xi[1,,k]   else xi[j,,k],
                                 if(nxy) xi[1,j,k,]  else if(ny) xi[1,,k,]  else xi[j,,k,],
                                 if(nxy) xi[1,j,k,,] else if(ny) xi[1,,k,,] else xi[j,,k,,],
                                 stop("wrong length(ovar)"))
                    ## panel
                    panel(z, y=y., col=v.col[k], ...)
                }
            } else { # no v (still n.sim or not)
                ## y values
                y. <- switch(length(d)-1, # length(ovar) - has.n.sim,
			     if(nxy) xi[1,j]   else if(ny) xi[1,]   else xi[j,],
			     if(nxy) xi[1,j,]  else if(ny) xi[1,,]  else xi[j,,],
			     if(nxy) xi[1,j,,] else if(ny) xi[1,,,] else xi[j,,,],
			     stop("wrong length(ovar)"))
                ## panel
		panel(z, y=y., col=v.col, ...)
            }

	    ## axes ############################################################

	    c.ax <- auxcol[1]

            ## x axis
	    if(i == n.y)
            {
		if(is.character(xv))
		    axis(1, at=z, labels=xv,
			 lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax)
		else
		    axis(1, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax)
                if(is.null(xlab)) xlab <- vList[[xvar]]$expr
		if(!identical(xlab, NA)) { # x axis label
		    ## would need "with local viewport": but the push..()/up..()
		    ## does not work  WHY ???
		    ## pushViewport(viewport(layout.pos.row=ny., layout.pos.col=j.))
		    ## grid.text(xlab, x=1, y=0.9, # check.overlap=TRUE,
		    ##           gp = gpar(col = c.ax))
		    ## popViewport()

		    ## Using traditional graphics instead:
		    ## mtext(xlab, side=1, at=zrange[2], col = c.ax, adj=c(0,0))
		    ## ... TODO MM ... use grid!  This *IS* ugly
		    ## draw it between the last two tick marks :
		    axt <- axTicks(1)
		    at <- mean(axt[-1:0 + length(axt)])
		    mtext(xlab, side=1, at=at, col = c.ax, font = 2)#, adj=0.5
		    ## draw axis arrow  '-->'  above the xlab:
		    ## pu <- par("usr"); if(par("xlog")) pu[1:2] <- 10^pu[1:2]
		    ## x0 <- pu[2] + (pu[2]-pu[1])
		    ## x1 <- x0 + (x0 - pu[1])* spc[1]/Pxsp/2
		    ## y0 <- pu[3]; if(par("ylog")) y0 <- 10^y0
		    ## arrows(x0=x0, y0=y0, x1=x1, length=1/8, col = c.ax, xpd=NA)
		}
	    }

            ## y axis
	    if(j == 1)
            {
		at <- axisTicks(par("usr")[3:4], log=par("ylog"))
		if(packageVersion("sfsmisc") >= "1.0-21")
		    ## allow for adjusting colors of small ticks
		    eaxis(2, at=at, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax, las = 1,
			  small.args=list(col=NA, col.ticks=c.ax, col.axis=c.ax))
		else
		    eaxis(2, at=at, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax, las = 1)
	    }
	    popViewport()

            ## panel labels ####################################################

	    ## column (col.vars) panel labels
	    if(i == 1 && !is.null(col.vars)) {
		pushViewport(viewport(layout.pos.row=1, layout.pos.col=j.))
		grid.rect(gp=gpar(col=NA, fill=auxcol[2]))
		grid.text(tryExpr(dn[[col.vars]][j]), x=0.5, y=0.5, gp=gpar(cex=tx.cex))
		popViewport()
	    }

	    ## row (row.vars) panel labels
	    if(j == n.x && !is.null(row.vars)) {
		pushViewport(viewport(layout.pos.row=i.,
                                      layout.pos.col=if(has.n.sim) nx.-1 else nx.))
		grid.rect(gp=gpar(col=NA, fill=auxcol[2]))
		grid.text(tryExpr(dn[[row.vars]][i]), x=0.5, y=0.5, gp=gpar(cex=tx.cex),
			  rot=-90)
		popViewport()
	    }

	} # for(j ..)
	if(verbose) cat("\n")

    } # for(i ..)


    ## additionals #############################################################

    ## ylab
    if(!identical(ylab, NA)) { # may be language where is.na(.) warns
	max.row <- if(lv) ny.-2 else ny.-1
	stopifnot(max.row >= 2) # defensive programming
        ## main
	pushViewport(viewport(layout.pos.row=2:max.row, layout.pos.col=1))
	grid.text(ylab, x=unit(0.6, "lines"), just=c("left", "centre"), rot=90)
	popViewport()
    }

    ## ylab2 (n.sim)
    if(has.n.sim && do.n.sim) {
	max.row <- if(lv) ny.-2 else ny.-1
	stopifnot(max.row >= 2) # defensive programming
        ## main
	pushViewport(viewport(layout.pos.row=2:max.row, layout.pos.col=nx.))
        ## note: vLnsim = vList[["n.sim"]] (which is defined if has.n.sim)
        lab <- substitute(lab == val, list(lab=vLnsim[["expr"]], val=vLnsim[["value"]]))
	grid.text(lab, y=unit(0, "npc"), # placing it at last/bottom panel(s)
                  just=c("left", "centre"), # given that point, use this (x, y) adjustment
                  rot=90, gp=gpar(cex=leg.cex))
	popViewport()
    }

    ## legend (for all variables appearing in a single panel)
    if(do.legend)
    {
        ## main
        wran <- if(has.n.sim) 2:(nx.-2) else 2:(nx.-1) # width range
	pushViewport(viewport(layout.pos.row=ny., layout.pos.col=wran))
        grid.draw(legend.grob)
        popViewport()
    }
    invisible(gl)
}
##' { end } mayplot



## ### pure grid (+lattice)

## if(FALSE) {##============ mayplotGrob() not yet ===================================

## ## From R-help (19 March 2010), by Baptiste Auguie
## latticeGrob <- function(x, ...) grob(.p. = x, cl="lattice", ...)
## ##                                   ---
## ##' when grid.draw() does its battery of things, drawDetails() is one of them:
## drawDetails.lattice <- function(x, recording=FALSE) {
##     lattice:::plot.trellis(x$.p., newpage=FALSE)
##     ##                       ---
## }

## mayplotGrob <- function(x, vList, row.vars, col.vars, xvar, log = "",
##                         panelGrob = NULL,
##                         spc = c(0.04/max(1,n.x-1), 0.04/max(1,n.y-1)),
##                         auxcol = c("gray40", "gray78", "gray90", "white"),
##                         pcol = c("blue", "red", "orange"),
##                         xlab = NULL, ylab = NULL,
##                         ylim = "global", type = "l",
##                         doBoxp = has.n.sim,
##                         verbose = getOption("verbose"),
##                         vp.ex = 1, draw=TRUE, ...)
## {
##     stopifnot(require("lattice"))

##     ## auxiliary function for building names
##     getName <- function(x) vapply(x, `[[`, "", "name")
##     ## basics
##     if(missing(row.vars)) row.vars <- NULL
##     if(missing(col.vars)) col.vars <- NULL
##     namD <- names( dn <- dimnames(x) )
##     if(!all(namD %in% names(vList)))
## 	stop("names(dimnames(.)) should all appear in names(vList)")
##     d <- dim(x)
##     if(is.null(names(d))) names(d) <- namD else stopifnot(names(d) == namD)
##     has.n.sim <- "n.sim" %in% namD
##     if(has.n.sim)	 {
## 	vLnsim <- vList[["n.sim"]]
## 	stopifnot(vLnsim[["type"]] == "N",
## 		  (n.sim <- d[["n.sim"]]) == vLnsim[["value"]])
##     }
##     else if(doBoxp) stop("no 'n.sim' -- doBoxp=TRUE  is not possible")
##     ## the variable displayed in one plot (with different colors):
##     ovar <- c(xvar, col.vars, row.vars, if(has.n.sim) "n.sim")
##     v <- setdiff(namD, ovar)

##     ## TODO MM: great if this worked with length(v) == 2 {"n" & "d"}, and maybe > 2
##     ## (e.g. *overplot*: one full plot for each "n" )
##     if(doBoxp && n.sim <= 4)
## 	warning("box plots with ",n.sim," values are rarely useful")
##     else if(!doBoxp && has.n.sim && n.sim > 1)
## 	stop("has n.sim (= ",n.sim,") > 1, but no boxplots.",
## 	     "\n *** error bars not yet implemented!")
##     if((lv <- length(v)) > 1)
## 	stop("length(v)=",lv," not yet implemented; too many variables (per plot panel)")
##     if(lv == 0) {
## 	if(verbose) cat("no 'inner plot' variable 'v'\n")
## 	v <- NULL
## 	v.col <- pcol[1]
##     } else { ## (lv >= 1)
## 	stopifnot(1 <= (nv <- length(dn[[v]])))
## 	## if(nv > length(pcol))# not an error; can still work!
## 	##     warning("less colors in 'pcol' than values in ",v,
## 	## 	    " (", length(pcol)," < ", nv,")")
## 	v.col <- colorRampPalette(pcol, space="Lab")(nv) # colors for v
##     }
##     ## permute to know the component indices (v and others may be NULL !):
##     p.x <- c(row.vars, col.vars, v, if(has.n.sim) "n.sim", xvar)
##     x <- aperm(x, perm = p.x)
##     if(verbose) {
## 	cat("xvar, col.vars, row.vars, v, n.sim : those *not* NULL are:\n")
## 	print(c(xvar=xvar, col.vars=col.vars,
##                 row.vars=row.vars, v=v, n.sim=if(has.n.sim) n.sim))
## 	if(is.null(panelGrob)) cat("doBoxp= ",doBoxp,"; ", sep="")
## 	stopifnot(identical(dim(x), unname(d[p.x])))
## 	cat("dim(x {after aperm()}):\n"); print(d[p.x])
##     }
##     if(is.null(xlab)) # default: the expression from varlist
## 	xlab <- vList[[xvar]]$expr
##     xv <- vList[[xvar]]$value
##     z <- if(is.character(xv) || is.recursive(xv)) seq_along(xv) else {
## 	z <- tryCatch(as.numeric(xv), error=function(e)e)
## 	if(inherits(z, "error")) seq_along(xv) else z
##     }
##     zrange <- range(z) # for forcing the same x axis limits per row
##     if(doBoxp) zrange <- extendrange(zrange)

##     if(is.null(panelGrob)) {
## 	markNA <- function(x.na, col) {
##             ## FIXME this is still base-graphics
## 	    xaxisGrob(1, at=x.na, labels=rep.int("N", length(x.na)),
## 		 mgp = c(1.5, 0.6, 0),
## 		 tick=FALSE, cex.axis=0.9, col.axis=adjustcolor(col, 0.8))
## 	    if(verbose)
## 		message("NA's at x = ", paste(format(x.na, collapse=", ")))
## 	}
## 	panelGrob <-
## 	    if(doBoxp)
## 		function(z, y, col, nv, ...) {
## 		    stopifnot(is.matrix(y), (p <- length(z)) == ncol(y))
##                     stop("doBoxp  not yet implemented here")
## 		    wex <- 0.5 * diff(range(z))/ p
## 		    boxplot.matrix(y, at = z, boxwex = wex, col=NA, border = col,
## 				   axes=FALSE, add=TRUE, names=FALSE, ...)
## 		    my <- apply(y, 2, median, na.rm=TRUE)
## 		    lines(z, my, col=col) # <- not yet '...' (warnings)
## 		    if(any(iy <- is.na(y))) { ## Mark  NA's at the bottom
## 			markNA(z[apply(iy, 2, any)], col)
## 		    }
## 		}
## 	    else
## 		function(z, y, col, nv, ...) {
## 		    ## plot corresponding points/lines
##                     lg <- linesGrob(z, y, type=type, col=col, ...)
## 		    if(any(iy <- is.na(y))) { ## Mark NA's at the bottom
## 			markNA(z[iy], col)
## 		    }
## 		}

##         ## Workaround for now
##         panelGrob <- function(z, y, col, nv, ...)
##             rectGrob(gp = gpar(fill = "tomato")) # TODO for now!

##     }
##     else stopifnot(is.function(panelGrob),
## 		     length(ff <- formals(panelGrob)) >= 4,
## 		     ## the first argument can be named arbitrarily; the others must match:
## 		     c("y", "col","nv") %in% names(ff[-1]))

##     ylim.msg <- "'ylim' must be \"global\", \"local\", or numeric of length two"
##     if(is.character(ylim)) {
## 	ylim <- match.arg(ylim, c("global", "local"))
## 	do.ylim <- ylim == "local"
## 	if(!do.ylim) ylim <- range(x, na.rm=TRUE)
##     } else if(is.numeric(ylim)) {
## 	do.ylim <- FALSE
## 	if(!(length(ylim) == 2 && is.finite(ylim)))
## 	    stop(ylim.msg)
##     } else stop(ylim.msg)

##     ## set up the grid layout
##     getdn <- function(nm) if(is.null(nm)) nm else dn[[nm]]
##     nx <- length(getdn(col.vars)) # number of plot columns
##     ny <- length(getdn(row.vars)) # number of plot rows

##     n.y <- max(1, ny)
##     n.x <- max(1, nx)
##     stopifnot(length(spc) == 2)

##     ## build matrix structure

##     ## basic matrix layout for panels
##     res <- gtable(widths =unit(rep(1, n.x), "null"),
##                   heights=unit(rep(1, n.y), "null"))
##     ## gmat <- gtable(## +1: add 1 column at the right border for 'row.var's:
##     ##                widths  = unit(rep(1, n.x + 1), "null"),
##     ##                ## +1: add 1 row    at the top   border for 'col.var's:
##     ##                heights = unit(rep(1, n.y + 1), "null"))

## ###=========== WORK FROM HERE ========================================

##     ## panels ##################################################################

##     npanel <- function(ji)
##     {
##         i <- ji[2]
## 	xi <-
## 	    if(ny) switch(length(d),
## 			  x[i ,drop=FALSE],
## 			  x[i , ,drop=FALSE],
## 			  x[i ,, ,drop=FALSE],
## 			  x[i ,,, ,drop=FALSE],
## 			  x[i ,,,, ,drop=FALSE],
## 			  stop(gettextf("'x' has wrong rank = length(dim(x)) = %d",
##                                         length(d))))
## 	    else x

##         j <- ji[1]
## 	if(do.ylim)
## 	    ylim <- range(xi, na.rm=TRUE) # for forcing the same y axis limits per row
## 	    ## j. <- 2*j # row index in layout (for jumping over gaps)
## 	    if(verbose) cat(sprintf("%d (%d) ", j))#, j.))

##         ## prelim stuff
##         ch <- paste0("panel (", io, ",", j, ")")
##         nxy <- nx && ny
##         if(lv) {
##             ## FIXME: for(k in 1:nv)
##             warning("nv=",nv, "; we are only drawing the first")
##             k <- 1
##             x. <- ## panel(z, y=
##                   switch(length(d)-2,   #length(ovar)-1 -has.n.sim,
##                          if(nxy) xi[1,j,k]   else if(ny) xi[1,,k]   else xi[j,,k],
##                          if(nxy) xi[1,j,k,]  else if(ny) xi[1,,k,]  else xi[j,,k,],
##                          if(nxy) xi[1,j,k,,] else if(ny) xi[1,,k,,] else xi[j,,k,,],
##                          stop("wrong length(ovar)"))
##                   ## ,
##                   ## col = v.col[k], nv=nv, ...)
##         } else # no 'v';  (still n.sim or not)
##             x. <- ## panel(z, y=
##                   switch(length(d)-1, #length(ovar) -has.n.sim,
##                          if(nxy) xi[1,j]   else if(ny) xi[1,]   else xi[j,],
##                          if(nxy) xi[1,j,]  else if(ny) xi[1,,]  else xi[j,,],
##                          if(nxy) xi[1,j,,] else if(ny) xi[1,,,] else xi[j,,,],
##                          stop("wrong length(ovar)"))
##                   ## ,
##                   ## col = v.col, nv=0, ...)

##         ## r.x <- length(dx <- dim(x))
##         ## if(r.x == 5)
##         ##     x. <- x[i,j,,1,1] # points (for fixed 'tau', 'family'; additionally, we fix 'n' and 'd' here)
##         ## else if(r.x == 4)
##         ##     x. <- x[i,j,1, ] # points ...
##         ## else if(r.x == 3)
##         ##     x. <- x[i,j, ] # points ...
##         ## else stop("silly panel not yet working for length(dim(x)) = ", r.x)
##         ## plot x. (in alpha for each d)
##         #z <- c(0.9, 0.95, 0.99, 0.999)
##         xyscales <- list(y=list(log=TRUE, equispaced.log = FALSE),
##                          draw = FALSE)
## 	xy <- xyplot(x. ~ z, ## 2^(11:30)~11:30,
##                      scales= xyscales)
##         gxy <- latticeGrob(xy, name="xy-grob")
##         gTree(rectGrob(gp=gpar(fill=auxcol[3])), # background

##                  ## __TODO__: plot.window(zrange, ylim, log=log)

##                  ## gridGrob(.........)
##                  ## grid(col=auxcol[4], lty="solid", lwd=grid.lwd, equilogs=FALSE)

##                  if(verbose) textGrob(label=paste("Hi, I'm", ch)),
##                  if(FALSE)
##                  linesGrob(x = z,#c(0.9, 0.95, 0.99, 0.999),
##                            y = x.,
##                            default.units = "native",
##                            vp = dataViewport(z, x.)),
##                  gxy,
##                  name=ch)
##     }

##     ## construct panels
##     ji.mat <- expand.grid(col=seq_len(n.x), row=seq_len(n.y),
##                           KEEP.OUT.ATTRS = FALSE)
##     panels <- apply(ji.mat, 1, npanel) # compute panels
##     ## add panels to plot (~> existing layout)
##     ## (t, l): = (top,left)-positions (usual matrix index) in the layout;
##     ## see res$layout and source of gtable_add_grob
##     res <- gtable_add_grob(res, grobs=panels,
##                            t = ji.mat[,2], # 1, ..., 1, 2, ..., 2, n.y, ..., n.y
##                            l = ji.mat[,1], # 1, 2, ..., n.x, 1, 2, ..., n.x
##                            name=getName(panels))


##     ## ## for y axis label; (nx-1)* (plot+gaps); plot; row.vars labels :
##     ## ## space for 1 (of n.x) plot in x direction:
##     ## Pxsp <- (1 - (axlabspc[1] + labspc[1]+ (n.x-1)*spc[1]))/ n.x
##     ## wid.gl <- c(axlabspc[1], rep(c(Pxsp, spc[1]), n.x-1), Pxsp, labspc[1])
##     ## ## for col.vars labels; (ny-1)*(plot+gap); plot; if(lv): {legend for v}
##     ## ## space for 1 (of n.y) plot in y direction:
##     ## leg.sp <- if(lv) legspc else 0
##     ## Pysp <- (1 - (labspc[2]+ axlabspc[2]+ leg.sp + (n.y-1)*spc[2]))/ n.y
##     ## hei.gl <- c(labspc[2], rep(c(Pysp, spc[2]), n.y-1), Pysp, axlabspc[2],
##     ##     	if(lv) leg.sp)
##     ## nx. <- length(wid.gl)
##     ## ny. <- length(hei.gl)

##     ## ## plot settings, restored on exit
##     ## opar <- par(no.readonly=TRUE); on.exit(par(opar))
##     ## plot.new() # start (empty) new page with 'graphics'
##     ## ## TODO MM: use smart units for widths/heights rather than just 'npc'
##     ## gl <- grid.layout(nrow=ny., ncol=nx., widths= wid.gl, heights = hei.gl,
##     ##          ## units in npc as for pdf(); no square plotting region otherwise:
##     ##          default.units="npc")

##     ## if(show.layout) grid.show.layout(gl, vp=viewport(width=1.25, height=1.25))
##     ## pushViewport(viewport(layout=gl)) # use this layout in a viewport


## ###---- currently  "muted":
## if(FALSE)
##     ## --- plot data ---
##     for(i in 1:n.y) { # rows
## 	## i. <- 2*i # column index in layout (for jumping over gaps)
## 	if(verbose) cat(sprintf("plot row %d (%d): [columns:] ", i))#, i.))
## 	xi <-
## 	    if(ny) switch(length(d),
## 			  x[i ,drop=FALSE],
## 			  x[i , ,drop=FALSE],
## 			  x[i ,, ,drop=FALSE],
## 			  x[i ,,, ,drop=FALSE],
## 			  x[i ,,,, ,drop=FALSE],
## 			  stop(gettextf("'x' has wrong rank = length(dim(x)) = %d",
##                                         length(d))))
## 	    else x
## 	if(do.ylim)
## 	    ylim <- range(xi, na.rm=TRUE) # for forcing the same y axis limits per row
## 	for(j in 1:n.x) { # columns
## 	    ## j. <- 2*j # row index in layout (for jumping over gaps)
## 	    if(verbose) cat(sprintf("%d (%d) ", j))#, j.))
## 	    ## pushViewport(viewport(layout.pos.row=i., layout.pos.col=j.))

## 	    ## plot
## 	    ## grid.rect(gp=gpar(col=NA, fill=auxcol[3])) # background
## 	    ## ## start a 'graphics' plot
## 	    ## par(plt = gridPLT())
## 	    ## par(new=TRUE) # always do this before each new 'graphics' plot
## 	    ## ## set up coordinate system:
## 	    ## plot.window(zrange, ylim, log=log)
## 	    ## ## background grid:
## 	    ## grid(col=auxcol[4], lty="solid", lwd=grid.lwd, equilogs=FALSE)

## 	    nxy <- nx && ny
## 	    if(lv) for(k in 1:nv)
## 		panel(z, y=
## 		      switch(length(d)-2, #length(ovar)-1 -has.n.sim,
## 			     if(nxy) xi[1,j,k]   else if(ny) xi[1,,k]   else xi[j,,k],
## 			     if(nxy) xi[1,j,k,]  else if(ny) xi[1,,k,]  else xi[j,,k,],
## 			     if(nxy) xi[1,j,k,,] else if(ny) xi[1,,k,,] else xi[j,,k,,],
## 			     stop("wrong length(ovar)"))
## 		      ,
## 		      col = v.col[k], nv=nv, ...)
## 	    else # no 'v';  (still n.sim or not)
## 		panel(z, y=
## 		      switch(length(d)-1, #length(ovar) -has.n.sim,
## 			     if(nxy) xi[1,j]   else if(ny) xi[1,]   else xi[j,],
## 			     if(nxy) xi[1,j,]  else if(ny) xi[1,,]  else xi[j,,],
## 			     if(nxy) xi[1,j,,] else if(ny) xi[1,,,] else xi[j,,,],
## 			     stop("wrong length(ovar)"))
## 		      ,
## 		      col = v.col, nv=0, ...)
## 	    ## axes
## 	    c.ax <- auxcol[1]
## 	    if(i == n.y) { # x axis
## 		if(is.character(xv))
## 		    axis(1, at=z, labels=xv,
## 			 lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax)
## 		else
## 		    axis(1, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax)
## 		if(!is.null(xlab) && !identical(xlab, NA)) { ## x axis label
## 		    ## would need  "with local viewport": but the push..()/up..()
## 		    ## does not work  WHY ???
## 		    ## pushViewport(viewport(layout.pos.row=ny., layout.pos.col=j.))
## 		    ## grid.text(xlab, x=1, y=0.9, # check.overlap=TRUE,
## 		    ##           gp = gpar(col = c.ax))
## 		    ## popViewport()

## 		    ## Using traditional graphics instead:
## 		    ## mtext(xlab, side=1, at=zrange[2], col = c.ax, adj=c(0,0))
## 		    ## ... TODO MM ... use grid!  This *IS* ugly
## 		    ## draw it between the last two tick marks :
## 		    axt <- axTicks(1)
## 		    at <- mean(axt[-1:0 + length(axt)])
## 		    mtext(xlab, side=1, at=at, col = c.ax, font = 2)#, adj=0.5
## 		    ## draw axis arrow  '-->'  above the xlab:
## 		    ## pu <- par("usr"); if(par("xlog")) pu[1:2] <- 10^pu[1:2]
## 		    ## x0 <- pu[2] + (pu[2]-pu[1])
## 		    ## x1 <- x0 + (x0 - pu[1])* spc[1]/Pxsp/2
## 		    ## y0 <- pu[3]; if(par("ylog")) y0 <- 10^y0
## 		    ## arrows(x0=x0, y0=y0, x1=x1, length=1/8, col = c.ax, xpd=NA)
## 		}
## 	    }
## 	    if(j == 1) { # y axis
## 		at <- axisTicks(par("usr")[3:4], log=par("ylog"))
## 		if(packageVersion("sfsmisc") >= "1.0-21")
## 		    ## allow for adjusting colors of small ticks
## 		    eaxis(2, at=at, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax,
## 			  las = ylas,
## 			  small.args=list(col=NA, col.ticks=c.ax, col.axis=c.ax))
## 		else
## 		    eaxis(2, at=at, lwd=ax.lwd, col=NA, col.ticks=c.ax, col.axis=c.ax,
## 			  las = ylas)
## 	    }
## 	    popViewport()

## 	    ## column labels (col.vars)
## 	    if(i == 1 && !is.null(col.vars)) {
## 		pushViewport(viewport(layout.pos.row=1, layout.pos.col=j.))
## 		grid.rect(gp=gpar(col=NA, fill=auxcol[2]))
## 		grid.text(tryExpr(dn[[col.vars]][j]), x=0.5, y=0.5, gp=gpar(cex=tx.cex))
## 		popViewport()
## 	    }

## 	    ## row labels (row.vars)
## 	    if(j == n.x && !is.null(row.vars)) {
## 		pushViewport(viewport(layout.pos.row=i., layout.pos.col=nx.))
## 		grid.rect(gp=gpar(col=NA, fill=auxcol[2]))
## 		grid.text(tryExpr(dn[[row.vars]][i]), x=0.5, y=0.5, gp=gpar(cex=tx.cex),
## 			  rot=-90)
## 		popViewport()
## 	    }
## 	}## for(j ..)
## 	if(verbose) cat("\n")
##     }## for(i ..)
## ###---- currently  "muted"

##     ## panel labels ############################################################

##     panelLabs <- function(i, rot=0, base.name=if(rot) "row" else "col")
##         gTree(rectGrob(gp=gpar(fill="gray50")),
##               textGrob(label=paste(base.name, "lab", i), rot=rot),
##               name=paste("label", base.name, i))

##     col.pLabs <- lapply(seq_len(n.x), panelLabs)          # list of col panel label grobs
##     row.pLabs <- lapply(seq_len(n.y), panelLabs, rot=-90) # list of row panel label grobs

##     ## add col panel labels
##     ## layout for col panel labels: pos = 0: on top
##     res <- gtable_add_rows(res, heights=2*grobHeight(textGrob("X")), pos=0)
##     res <- gtable_add_grob(res, grobs=col.pLabs, t=1, l=seq_len(n.x),
##                            name=getName(col.pLabs)) # fill in labels
##     ## add row panel labels
##     ## layout for row panel labels (pos = default): at "end", i.e., bottom:
##     res <- gtable_add_cols(res, widths=2*grobHeight(textGrob("X")))
##     res <- gtable_add_grob(res, grobs=row.pLabs, t=1+seq_len(n.y), l=n.x+1,
##                            name=getName(row.pLabs)) # fill in labels
##     ## note: the t=1+... is due to the fact that we already added a row on top of the basic layout

##     ## inter-row and -column space #############################################

##     as.unit <- function(x) if(is.unit(x)) x else unit(x, "npc")
##     spc <- as.unit(spc)
##     for(i in rev(1L +seq_len(n.y-1L)))# "1+" already have top "col lab" labels
## 	res <- gtable_add_rows(res, heights= spc[2], pos=i) # row space
##     for(j in rev(seq_len(n.x-1L)))
## 	res <- gtable_add_cols(res, widths = spc[1], pos=j) # col space


##     ## labels ##################################################################

##     res <- gtable_add_cols(res, widths =2*grobHeight(textGrob("X")))
##     res <- gtable_add_rows(res, heights=2*grobHeight(textGrob("X")))

##     ## TODO continue... (with filling in labels)

##     if(draw) grid.draw(res)
##     ## return
##     invisible(res)
## }

## }##============ mayplotGrob() not yet ===================================
