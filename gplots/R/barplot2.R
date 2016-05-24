# Revision 2.1 2005/06/06
# - Modified default behavior with 0's and NA's in
#   'height' so that these values are not plotted.
# - Warning messages added in the case of the above.

# Revision 2.0 2005/04/27
# - Added panel.first and panel.last arguments
# - As per R 2.0.0, the default barplot() method by default uses a
#   gamma-corrected grey palette (rather than the heat color
#   palette) for coloring its output when given a matrix.

barplot2 <- function(height, ...) UseMethod("barplot2")

barplot2.default <-
    function(
        height,
        width = 1,
        space = NULL,
        names.arg = NULL,
        legend.text = NULL,
        beside = FALSE,
        horiz = FALSE,
        density = NULL,
        angle = 45,
        col = NULL,
        prcol = NULL,
        border = par("fg"),
        main = NULL,
        sub = NULL,
         xlab = NULL,
         ylab = NULL,

        xlim = NULL,
        ylim = NULL,
        xpd = TRUE,
        log = "",

        axes = TRUE,
        axisnames = TRUE,

        cex.axis = par("cex.axis"),
        cex.names = par("cex.axis"),

        inside = TRUE,
        plot = TRUE,
        axis.lty = 0,
        offset = 0,

        plot.ci = FALSE,
        ci.l = NULL,
        ci.u = NULL,

        ci.color = "black",
        ci.lty = "solid",
        ci.lwd = 1,
        ci.width = 0.5,

        plot.grid = FALSE,
        grid.inc = NULL,

        grid.lty = "dotted",
        grid.lwd = 1,
        grid.col = "black",

        add = FALSE,
        panel.first = NULL,
        panel.last = NULL,
        ...)
{
    if (!missing(inside)) .NotYetUsed("inside", error = FALSE)# -> help(.)

    if (missing(space))
      space <- if (is.matrix(height) && beside) c(0, 1) else 0.2
    space <- space * mean(width)

    if (plot && axisnames && missing(names.arg))
      names.arg <-
          if(is.matrix(height)) colnames(height) else names(height)

    if (is.vector(height)
        || (is.array(height) && (length(dim(height)) == 1))) {
        ## Treat vectors and 1-d arrays the same.
	height <- cbind(height)
	beside <- TRUE
        ## The above may look strange, but in particular makes color
        ## specs work as most likely expected by the users.
        if(is.null(col)) col <- "grey"
    } else if (is.matrix(height)) {
        ## In the matrix case, we use " heat colors" by default.
        if(is.null(col)) col <- heat.colors(nrow(height))
    }
    else
	stop(paste(sQuote("height"), "must be a vector or a matrix"))

    if(is.logical(legend.text))
      legend.text <-
        if(legend.text && is.matrix(height)) rownames(height)

    # Check for log scales
    logx <- FALSE
    logy <- FALSE

    if (log != "")
    {
      if (any(grep("x", log)))
        logx <- TRUE

      if (any(grep("y", log)))
        logy <- TRUE
    }

    # Cannot "hatch" with rect() when log scales used
    if ((logx || logy) && !is.null(density))
      stop("Cannot use shading lines in bars when log scale is used")

    NR <- nrow(height)
    NC <- ncol(height)

    if (beside) {
	if (length(space) == 2)
	    space <- rep.int(c(space[2], rep.int(space[1], NR - 1)), NC)
	width <- rep(width, length.out = NR)
    } else
	width <- rep(width, length.out = NC)

    offset <- rep(as.vector(offset), length.out = length(width))

    delta <- width / 2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta

    #if graphic will be stacked bars, do not plot ci
    if (!beside && (NR > 1) && plot.ci)
      plot.ci = FALSE

    # error check ci arguments
    if (plot && plot.ci)
    {
      if ((missing(ci.l)) || (missing(ci.u)))
        stop("confidence interval values are missing")

      if (is.vector(ci.l)
        || (is.array(ci.l) && (length(dim(ci.l)) == 1)))
        ci.l <- cbind(ci.l)
      else if (!is.matrix(ci.l))
        stop(paste(sQuote("ci.l"), "must be a vector or a matrix"))

      if (is.vector(ci.u)
        || (is.array(ci.u) && (length(dim(ci.u)) == 1)))
        ci.u <- cbind(ci.u)
      else if (!is.matrix(ci.u))
        stop(paste(sQuote("ci.u"), "must be a vector or a matrix"))

      if (any(dim(height) != dim(ci.u)))
        stop(paste(sQuote("height"), "and", sQuote("ci.u"),
                   "must have the same dimensions."))
      else if (any(dim(height) != dim(ci.l)))
        stop(paste(sQuote("height"), "and", sQuote("ci.l"),
                   "must have the same dimensions."))
    }

    # check height + offset/ci.l if using log scale to prevent log(<=0) error
    # adjust appropriate ranges and bar base values
    if ((logx && horiz) || (logy && !horiz))
    {
      # Check for NA values and issue warning if required
      height.na <- sum(is.na(height))
      if (height.na > 0)
      {
         warning(sprintf("%.0f values == NA in 'height' omitted from logarithmic plot",
                          height.na), domain = NA)
      }

      # Check for 0 values and issue warning if required
      # _FOR NOW_ change 0's to NA's so that other calculations are not
      # affected. 0's and NA's affect plot output in the same way anyway,
      # except for stacked bars, so don't change those.
      height.lte0 <- sum(height <= 0, na.rm = TRUE)
      if (height.lte0 > 0)
      {
        warning(sprintf("%0.f values <=0 in 'height' omitted from logarithmic plot",
                         height.lte0), domain = NA)

        # If NOT stacked bars, modify 'height'
        if (beside)
          height[height <= 0] <- NA
      }

      if (plot.ci && (min(ci.l) <= 0))
        stop("log scale error: at least one lower c.i. value <= 0")

      if (logx && !is.null(xlim) && (xlim[1] <= 0))
        stop("log scale error: 'xlim[1]' <= 0")

      if (logy && !is.null(ylim) && (ylim[1] <= 0))
        stop("'log scale error: 'ylim[1]' <= 0")

      # arbitrary adjustment to display some of bar for min(height) since
      # 0 cannot be used with log scales. If plot.ci, also check ci.l
      if (plot.ci)
      {
        rectbase <- c(height[is.finite(height)], ci.l)
        rectbase <- min(0.9 * rectbase[rectbase > 0])
      }
      else
      {
        rectbase <- height[is.finite(height)]
        rectbase <- min(0.9 * rectbase[rectbase > 0])
      }

      # if axis limit is set to < above, adjust bar base value
      # to draw a full bar
      if (logy && !is.null(ylim) && !horiz)
        rectbase <- ylim[1]
      else if (logx && !is.null(xlim) && horiz)
        rectbase <- xlim[1]

      # if stacked bar, set up base/cumsum levels, adjusting for log scale
      if (!beside)
        height <- rbind(rectbase, apply(height, 2, cumsum))

      # if plot.ci, be sure that appropriate axis limits are set to include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height

      rangeadj <- c(0.9 * lim + offset, lim + offset)
      rangeadj <- rangeadj[rangeadj > 0]
    }
    else
    {
      # Use original bar base value
      rectbase <- 0

      # if stacked bar, set up base/cumsum levels
      if (!beside)
        height <- rbind(rectbase, apply(height, 2, cumsum))

      # if plot.ci, be sure that appropriate axis limits are set to include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height

      # use original range adjustment factor
      rangeadj <- c(-0.01 * lim + offset, lim + offset)
    }

    # define xlim and ylim, adjusting for log-scale if needed
    if (horiz)
    {
      if (missing(xlim)) xlim <- range(rangeadj, na.rm=TRUE)
      if (missing(ylim)) ylim <- c(min(w.l), max(w.r))
    }
    else
    {
      if (missing(xlim)) xlim <- c(min(w.l), max(w.r))
      if (missing(ylim)) ylim <- range(rangeadj, na.rm=TRUE)
    }

    if (beside)
      w.m <- matrix(w.m, ncol = NC)

    if(plot) ##-------- Plotting :
    {
      opar <-
        if (horiz) par(xaxs = "i", xpd = xpd)
        else       par(yaxs = "i", xpd = xpd)
      on.exit(par(opar))

      # If add = FALSE open new plot window
      # else allow for adding new plot to existing window
      if (!add)
      {
        plot.new()
        plot.window(xlim, ylim, log = log, ...)
      }

      # Execute the panel.first expression. This will work here
      # even if 'add = TRUE'
      panel.first

      # Set plot region coordinates
      usr <- par("usr")

      # adjust par("usr") values if log scale(s) used
      if (logx)
      {
        usr[1] <- 10 ^ usr[1]
        usr[2] <- 10 ^ usr[2]
      }

      if (logy)
      {
        usr[3] <- 10 ^ usr[3]
        usr[4] <- 10 ^ usr[4]
      }

      # if prcol specified, set plot region color
      if (!missing(prcol))
        rect(usr[1], usr[3], usr[2], usr[4], col = prcol)

      # if plot.grid, draw major y-axis lines if vertical or x axis if horizontal
      # R V1.6.0 provided axTicks() as an R equivalent of the C code for
      # CreateAtVector.  Use this to determine default axis tick marks when log
      # scale used to be consistent when no grid is plotted.
      # Otherwise if grid.inc is specified, use pretty()

      if (plot.grid)
      {
        par(xpd = FALSE)

        if (is.null(grid.inc))
        {
          if (horiz)
          {
            grid <- axTicks(1)
            abline(v = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
          }
          else
          {
            grid <- axTicks(2)
            abline(h = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
          }
        }
        else
        {
          if (horiz)
          {
            grid <- pretty(xlim, n = grid.inc)
            abline(v = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
          }
          else
          {
            grid <- pretty(ylim, n = grid.inc)
            abline(h = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
          }
        }

         par(xpd = xpd)
      }

      xyrect <- function(x1,y1, x2,y2, horizontal = TRUE, ...)
      {
        if(horizontal)
          rect(x1,y1, x2,y2, ...)
        else
          rect(y1,x1, y2,x2, ...)
      }

      if (beside)
        xyrect(rectbase + offset, w.l, c(height) + offset, w.r, horizontal=horiz,
               angle = angle, density = density, col = col, border = border)
      else
      {
        for (i in 1:NC)
          xyrect(height[1:NR, i] + offset[i], w.l[i], height[-1, i] + offset[i], w.r[i],
                 horizontal=horiz, angle = angle, density = density,
                 col = col, border = border)
      }

      # Execute the panel.last expression here
      panel.last

      if (plot.ci)
      {
        # CI plot width = barwidth / 2
        half.ci.width = width * ci.width / 2

        if (horiz)
        {
          segments(ci.l, w.m, ci.u, w.m, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(ci.l, w.m - half.ci.width, ci.l, w.m + half.ci.width, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(ci.u, w.m - half.ci.width, ci.u, w.m + half.ci.width, col = ci.color, lty = ci.lty, lwd = ci.lwd)
        }
        else
        {
          segments(w.m, ci.l, w.m, ci.u, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(w.m - half.ci.width, ci.l, w.m + half.ci.width, ci.l, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(w.m - half.ci.width, ci.u, w.m + half.ci.width, ci.u, col = ci.color, lty = ci.lty, lwd = ci.lwd)
        }
      }

      if (axisnames && !is.null(names.arg)) # specified or from {col}names
      {
        at.l <-
        if (length(names.arg) != length(w.m))
        {
          if (length(names.arg) == NC) # i.e. beside (!)
            colMeans(w.m)
          else
            stop("incorrect number of names")
        }
        else w.m

        axis(if(horiz) 2 else 1, at = at.l, labels = names.arg, lty = axis.lty, cex.axis = cex.names, ...)
      }

      if(!is.null(legend.text))
      {
        legend.col <- rep(col, length = length(legend.text))

        if((horiz & beside) || (!horiz & !beside))
        {
          legend.text <- rev(legend.text)
          legend.col <- rev(legend.col)
          density <- rev(density)
          angle <- rev(angle)
        }

        # adjust legend x and y values if log scaling in use
        if (logx)
          legx <- usr[2] - ((usr[2] - usr[1]) / 10)
        else
          legx <- usr[2] - xinch(0.1)

        if (logy)
          legy <- usr[4] - ((usr[4] - usr[3]) / 10)
        else
          legy <- usr[4] - yinch(0.1)

        legend(legx, legy,
        legend = legend.text, angle = angle, density = density,
        fill = legend.col, xjust = 1, yjust = 1)
      }

      title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)

      # if axis is to be plotted, adjust for grid "at" values
      if (axes)
      {
        if(plot.grid)
          axis(if(horiz) 1 else 2, at = grid, cex.axis = cex.axis, ...)
        else
          axis(if(horiz) 1 else 2, cex.axis = cex.axis, ...)
      }

      invisible(w.m)

    }

    else w.m
}
