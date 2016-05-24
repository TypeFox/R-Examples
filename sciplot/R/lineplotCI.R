lineplot.CI <-
  function(x.factor, response, group = NULL, type = "b", x.cont=FALSE,
           legend = TRUE, trace.label = NULL, leg.lab = NULL,
           fixed = FALSE, x.leg = NULL, y.leg = NULL, cex.leg = 1, ncol = 1,
           pch = c(16,21,15,22,17,24, c(3:14)),
           fun = function(x) mean(x, na.rm=TRUE),
           ci.fun = function(x) c(fun(x)-se(x), fun(x)+se(x)),
           err.width = if(length(levels(as.factor(x.factor)))>10) 0 else .1,
           err.col = col, err.lty = 1, xlim = NULL, ylim = NULL,
           cex = NULL, lwd = NULL, col = "black", cex.axis = 1,
           xaxt = "s", data = NULL, subset = NULL, ...) {
  # Set up environment
  subset <- eval(substitute(subset), envir=data)
  fun <- eval(substitute(fun), envir=data)
  if(!is.null(data)) {
    if(!is.null(subset)) data <- subset(data,subset)
    x.factor <- eval(substitute(x.factor), envir=data)
    response <- eval(substitute(response), envir=data)
    group <- eval(substitute(group), envir=data)
  }

  subset = NULL
  
  #////////////////////////////////////////////////////////////////////#
  # Below modified from the interaction.plot function in stats package #
  # Take the legend section out to allow greater flexibility           #
  #////////////////////////////////////////////////////////////////////#

  int.plot <- function (x.factor = x.factor, group = group,
                        response = response, type = c("l","p", "b"),
                        legend = legend,
                        trace.label = deparse(substitute(group)),
                        fixed = FALSE, xlab = deparse(substitute(x.factor)),
                        ylab = ylabel, lty = nc:1, pch = NA, xpd = NULL,
                        leg.bg = par("bg"), leg.bty = "n", xtick = FALSE,
                        xlim=xlim, ylim=ylim, axes = TRUE, ...) {
  
    ylabel <- paste(deparse(substitute(fun)), "of ", 
                    deparse(substitute(response)))
    type <- match.arg(type)
    cells <- tapply(response, list(x.factor, group), fun)
    nr <- nrow(cells)
    nc <- ncol(cells)
    xvals <- if(x.cont) as.numeric(levels(as.factor(x.factor))) else 1:nr
    if (is.ordered(x.factor)) {
      wn <- getOption("warn")
      options(warn = -1)
      xnm <- as.numeric(levels(x.factor))
      options(warn = wn)
      if (!any(is.na(xnm)))
        xvals <- xnm
    }
    xlabs <- rownames(cells)
    ylabs <- colnames(cells)
    nch <- max(sapply(ylabs, nchar, type = "width"))
    if (is.null(xlabs))
      xlabs <- as.character(xvals)
    if (is.null(ylabs))
      ylabs <- as.character(1:nc)
    if(is.null(xlim)) { 
      xlim <- range(xvals)
      xleg <- xlim[2] + 0.05 * diff(xlim)
      xlim <- xlim + c(-0.2/nr, if (legend&is.null(x.leg)&is.null(y.leg))
                       0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)}
    else {
      xlim
      xleg <- xlim[2] - 0.25 * diff(xlim)}
    matplot(xvals, cells, ..., type = type, xlim = xlim, ylim = ylim,
            xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col,
            lty = lty, lwd=lwd, pch = pch)
    if (axes && xaxt != "n") {
      axisInt <- function(x, main, sub, lwd, bg, log, asp,
                          ...) axis(1, x, ...)
      mgp. <- par("mgp")
      if (!xtick)
        mgp.[2] <- 0
      axisInt(1, at = xvals, labels = xlabs, tick = xtick,
              mgp = mgp., xaxt = xaxt, ...)
    }
    ord <- sort.list(cells[nr, ], decreasing = TRUE)
    if (legend) {
      yrng <- diff(ylim)
      yleg <- ylim[2] - 0.1 * yrng
      if (!is.null(xpd) || {
        xpd. <- par("xpd")
        !is.na(xpd.) && !xpd. && (xpd <- TRUE)
      }) {
        op <- par(xpd = xpd)
        on.exit(par(op))
      }
      text(xleg, ylim[2] - 0.05 * yrng, paste("  ", trace.label),
           adj = 0)
      if (!fixed) {
        ylabs <- ylabs[ord]
        lty <- lty[1 + (ord - 1)%%length(lty)]
        col <- col[1 + (ord - 1)%%length(col)]
        pch <- pch[ord]
      }
    }
    invisible()

    return.data<-if(legend)
      list(pch=pch,ord=ord,xleg=xleg,yleg=yleg,ylabs=ylabs,lty=lty,
           leg.bty=leg.bty,leg.bg=leg.bg,ord=ord,xvals=xvals,cells=cells)
    else list(pch=pch,ord=ord,xvals=xvals)
       return(return.data)
  }
  #////////////////////////////////////////////////////////////////////#
  #                         End Section                                #
  #////////////////////////////////////////////////////////////////////#

  # Figure out if we're dealing with 1 or more total groups (i.e., we always
  # have 1 "x.factor" and we may have 1 or more "group".
  if(is.null(group)) groups = factor(x.factor) else {   
    # If more than 1 "y-factor", combine for plotting
    if(length(group[[1]]) > 1)
      group <- factor(interaction(group, lex.order=TRUE))
    # "groups" defines the combination of "x.factor" and "group"
    group <- factor(group) # Do this to drop unused levels
    groups <- list(x.factor, group)
  }

  # Calculate mean and SE's
  mn.data <- tapply(response, groups, fun)
  CI.data <- tapply(response, groups, ci.fun)
  # Determine y-axis plot region
  plot.limits = c(
    min(c(unlist(mn.data), unlist(CI.data)), na.rm=TRUE),
    max(c(unlist(mn.data), unlist(CI.data)), na.rm=TRUE))

  # Draw lines
  if(is.null(group)) {
    nlevels.x <-
      if(x.cont) as.numeric(levels(as.factor(x.factor)))
      else 1:nrow(mn.data)
    plot(nlevels.x, mn.data, xaxt="n", type=type, col=col, pch=NA, cex=cex,
         cex.axis=cex.axis, lwd=lwd,
         xlim=if(is.null(xlim)) {
           c(min(nlevels.x)-.2,max(nlevels.x)+.2)
         }
         else xlim,
         ylim=if(is.null(ylim)) plot.limits else ylim, ...)
    if(xaxt!="n") axis(1,labels=names(mn.data),
         at=nlevels.x, cex.axis=cex.axis, ...)
  }
  else leg.vals <-
    int.plot(x.factor, group, response, type = type,
             xlim = xlim, ylim = if(is.null(ylim)) plot.limits else ylim,
             cex.axis = cex.axis, trace.label = trace.label, pch = NA,
             legend = legend, ...)
  # Draw CI's
  if(is.null(group)) {
    nlevels.x <-
      if(x.cont) as.numeric(levels(as.factor(x.factor)))
      else 1:nrow(mn.data)
    CI.seln <- !is.na(mn.data)
    CI.plot <- matrix(unlist(CI.data[CI.seln]),
                      nrow=sum(CI.seln), byrow=TRUE)
    arrows(nlevels.x[CI.seln],
           CI.plot[,1],
           nlevels.x[CI.seln],
           CI.plot[,2],
           angle = 90, col = err.col, length = err.width, code = 3,
           lwd = lwd, lty=err.lty)
  }
  else {
    nlevels.y <- ncol(mn.data)
    for(i in 1:nlevels.y) {
      CI.seln <- !is.na(mn.data)[,i]
      CI.plot <- matrix(unlist(CI.data[CI.seln,i]),
                        nrow=sum(CI.seln), byrow=TRUE)
      arrows(leg.vals$xvals[CI.seln],
             CI.plot[,1],
             leg.vals$xvals[CI.seln],
             CI.plot[,2],
             angle=90, length=err.width,
             col=if(length(err.col)>1) err.col[i] else err.col,
             lty=if(length(err.lty)>1) err.lty[i] else err.lty,
             code=3,lwd=lwd)
    }
  }

  # Draw points (Note: adding points at this point allows points to be in
  # the foreground)
  if(type %in% c("p", "b")) {
    if(is.null(group)) {
      nlevels.x <-
        if(x.cont) as.numeric(levels(as.factor(x.factor)))
        else 1:nrow(mn.data)
      points(nlevels.x, mn.data,pch=pch[1],bg="white",cex=cex,col=col)
    }
    else {
      nlevels.y<-dim(mn.data)[2]
      for(i in 1:nlevels.y)
        points(leg.vals$xvals, mn.data[,i],pch=pch[i],bg="white",
               col=if(length(col)>1) col[i] else col,cex=cex)
    }
  }

  #////////////////////////////////////////////////////////////////////#
  # Below modified from the interaction.plot function in stats package #
  #////////////////////////////////////////////////////////////////////#

    if(legend & !is.null(group)) {
    legend(x=if(is.null(x.leg)) leg.vals$xleg else x.leg,
           y=if(is.null(y.leg)) leg.vals$yleg else y.leg,
           legend = if(!is.null(leg.lab)) leg.lab else {
             if(fixed) levels(as.factor(unlist(group))) else leg.vals$ylabs
           },
           pch = if(type %in% c("p", "b")) {
             if(!fixed) pch[leg.vals$ord]
             else pch
             },
           col = if(type %in% c("p", "b")) {
             if(!fixed  & length(col)>1) col[leg.vals$ord]
             else col
             },
           lty = if (type %in% c("l", "b")) {
             if(fixed) leg.vals$lty[order(leg.vals$ord)]
             else leg.vals$lty
             },
           ncol=ncol, bty = leg.vals$leg.bty, bg = leg.vals$leg.bg,
           cex=cex.leg)
  }

  #////////////////////////////////////////////////////////////////////#
  #                      End Section                                   #
  #////////////////////////////////////////////////////////////////////#

  invisible(list(vals=mn.data, CI=CI.data))
}
