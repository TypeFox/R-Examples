#' Plot of balance across multiple strata
#'
#' The plot allows a quick visual comparison of the effect of different
#' stratification designs on the comparability of different
#' variables. This is not a replacement for the omnibus statistical test
#' reported as part of \code{\link{print.xbal}}. This plot does allow the
#' analyst an easy way to identify variables that might be the primary culprits
#' of overall imbalances and/or a way to assess whether certain important
#' covariates might be imbalanced even if the omnibus test reports that
#' the stratification overall produces balance.
#'
#' By default all variables and all strata are plotted. The scope
#' of the plot can be reduced by using the \code{\link{subset.xbal}} function to
#' make a smaller \code{xbal} object with only the desired variables or
#' strata.
#'
#' \code{\link{xBalance}} can produce several different summary statistics for
#' each variable, any of which can serve as the data for this plot. By default,
#' the standardized differences between treated and control units makes a good
#' choice as all variables are on the same scale. Other statistics can be
#' selected using the \code{statistic} argument.
#'
#' @param x An object returned by \code{\link{xBalance}}
#' @param xlab The label for the x-axis of the plot
#' @param statistic The statistic to plot. The default choice of standardized
#' difference is a good choice as it will have roughly the same scale for all
#' plotted variables.
#' @param absolute Convert the results to the absolute value of the statistic.
#' @param strata.labels A named vector of the from \code{c(strata1 = "Strata Label 1", ...)}
#' that maps the stratification schemes to textual labels.
#' @param variable.labels A named vector of the from \code{c(var1 = "Var Label1", ...)}
#' that maps the variables to textual labels.
#' @param groups A vector of group names for each variable in
#' \code{x$results}. By default, factor level variables will be
#' grouped.
#' @param ... additional arguments to pass to \code{\link{balanceplot}}
#' @seealso \code{\link{xBalance}}, \code{\link{subset.xbal}}, \code{\link{balanceplot}}
#' @example inst/examples/plot.xbal.R
#' @import abind
#' @export
plot.xbal <- function(x,
                      xlab = "Standardized Differences",
                      statistic = "std.diff",
                      # not implemented: thecols = rainbow(length(which.strata)),
                      # not implemented: thesymbols = c(19,22,23,24,25)[1:length(which.strata)],
                      absolute = FALSE,
                      strata.labels = NULL,
                      variable.labels = NULL,
                      groups = NULL,
                      ...) {

  x <- prepareXbalForPlot(x, statistic, absolute, strata.labels, variable.labels)

  if (is.null(groups)) {
    groups <- attr(x, "groups")
  }

  return(balanceplot(x, xlab = xlab, groups = groups, ...))

  ### NOT RUN: (but saving while we transition to the more general balanceplot function

  # nvars <- dim(theresults)[1]
  # nstrata <- dim(theresults)[2]
  # varlabels <- rownames(theresults)

  # ypos <- seq(length.out = nvars)
  # xrange <- range(theresults, na.rm = TRUE)
  # xrange <- xrange + xrange * adjustxaxis

  # ##Setup the margin, adjust for the lengths of the labels.
  # par(mar=mar,tck=tck,mgp=mgp) ##set default margins
  # mymai<-par('mai')

  # # when using the SVG device, strwidth throws fits, so we hope that the default mai[2] is good enough
  # if (names(dev.cur()) != "svg") {
  #   mymai[2]<-max(c(strwidth(varlabels ,units="inches"),mymai[2]))
  # }

  # ##Setup the plotting region
  # par(mai=mymai)

  # if(length(thecols)!=length(which.strata)){
  #   if(length(thecols)==1){
  #     thecols<-rep(thecols,length(which.strata))
  #   }
  #   if(length(thecols)>1){
  #     cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
  #   }
  # }
  # names(thecols)<-which.strata

  # if(length(thesymbols)!=length(which.strata)){
  #   if(length(thesymbols)==1){
  #     thesymbols<-rep(thesymbols,length(which.strata))
  #   }
  #   if(length(thesymbols)>1){
  #     cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
  #   }
  # }
  # names(thesymbols)<-which.strata
  #
  #
  # plot(xrange,range(ypos),axes=FALSE,pch=19,col="blue",
  #      ylab="",xlab=thexlab,type="n",...)
  # for(i in which.strata){
  #   points(theresults[,i],ypos,col=thecols[i],pch=thesymbols[i])
  #   }
  # if(segments&length(which.strata)>1){ ##segments are mainly useful for drawing the eye to changes in balance along a single variable across more than 1 stratification
  #   for(j in ypos){
  #     segments(min(theresults[j,]),j,
  #              max(theresults[j,]),j,col=gray(.7),lwd=.5)
  #   }
  # }
  # axis(1,at=pretty(seq(xrange[1],xrange[2],length=5)))
  # axis(2,labels=varlabels,at=ypos,las=2,tick=FALSE)
  # lines(c(0,0),range(ypos)+c(-.025*length(ypos),.025*length(ypos)),col="grey",lwd=1)
  # ##segments(0,min(ypos),0,max(ypos),col="grey",lwd=1)
  # if(legend){
  #   legend(x="topright",#xrange[1],ypos[ypos==max(ypos)],
  #          legend=thestratalabs,
  #          col=thecols,
  #          pch=thesymbols,
  #          bty="n")
  # }
}

# Internal function for turning an xBalance object into something for `balanceplot`
prepareXbalForPlot <- function(x,
                               statistic = "std.diff",
                               absolute = FALSE,
                               strata.labels = NULL,
                               variable.labels = NULL,
                               ...) {

  if (dim(x$results)[2] > 1) {
    # this means that the user is passing an xBalance object with more than one statistic
    # so we need to trim it down

    # but first we need to make sure the statistic exists
    if (!(statistic %in% dimnames(x$results)[[2]])) {
      stop("Unknown statistic: ", statistic)
    }
    x <- subset(x, stats = statistic)
  }

  origs <- attr(x$results, "originals")

  x <- adrop(x$results, drop = 2)


  if (!is.null(variable.labels)) {
    if (is.null(names(variable.labels))) {
      stop("Variable labels must be a named vector of the form c('var1' = 'Var One', ...)")
    }
    rownames(x) <- variable.labels[rownames(x)]
  }

  if (!is.null(strata.labels)) {
    if (is.null(names(strata.labels))) {
      stop("Strata labels must be a named vector of the form c('var1' = 'Var One', ...)")
    }
    colnames(x) <- strata.labels[colnames(x)]
  }

  if (absolute) {
    x <- abs(x)
  }

  mgrps <- origs %in% names(which(table(origs) > 1))
  origs[!mgrps] <- NA
  attr(x, "groups") <- origs

  return(x)
}

#' Create a plot of the balance on variables across different stratifications.
#'
#' This plotting function summarizes variable by stratification matrices. For
#' each variable (a row in the \code{x} argument), the values are under each
#' stratification (the columns of \code{x}) plotted on the same line.
#'
#' It is conventional to standardize the differences to common scale
#' (e.g.  z-scores), but this is not required. When \code{ordered} is
#' set to true, plotting will automatically order the data from
#' largest imbalance to smallest based on the first column of
#' \code{x}.
#'
#' You can fine tune the colors and shapes with the like named
#' arguments. Any other arguments to the \code{\link{points}} function
#' can be passed in a list as \code{points.args}. Likewise, you can
#' fine tune the segments between points with \code{segments.args}.
#'
#' @param x A matrix of variables (rows) by strata (columns).
#' @param ordered Should the variables be ordered from
#' most to least imbalance on the first statistic?
#' @param segments Should lines be drawn between points for each
#' variable?
#' @param colors Either a vector or a matrix of shape indicators
#' suitable to use as a \code{col} argument to the
#' \code{\link{points}} function. If the argument is a vector, the
#' length should be the same as the number of columns in \code{x}. If
#' the argument is a matrix, it should have the same dims as \code{x}.
#' @param shapes Either a vector or a matrix of shape indicators
#' suitable to use as a \code{pch} argument to the
#' \code{\link{points}} function. If the argument is a vector, the
#' length should be the same as the number of columns in \code{x}. If
#' the argument is a matrix, it should have the same dims as
#' \code{x}. The suggested vector has been selected to work with
#' RSVGTipsDevice tool tips.
#' @param segments.args A list of arguments to pass to the
#' \code{\link{segments}} function.
#' @param points.args A list of arguments to pass to the \code{\link{points}} function.
#' @param xlab The label of the x-axis of the plot.
#' @param xrange The range of x-axis. By default, it is 1.25 times the range of \code{x}.
#' @param groups A factor that indicates the group of each row in
#' \code{x}. Groups are printed under a common header.
#' @param tiptext If you are using the \code{RSVGTipsDevice} library for
#' rendering, you can include an array of the dimensions of x
#' with another dimension of length 2. For example, if there are 4
#' observations and 2 strata, the array should be 4 by 2 by 2. The
#' \code{tiptext[i, j, 1]} entry will be the first line of the tool
#' tip for the data in \code{x[i, j]}. Likewise for the second row of
#' the tool tip.
#' @param include.legend Should a legend be included?
#' @param legend.title An optional title to attach to the legend.
#' @param ... Additional arguments to pass to \code{\link{plot.default}}.
#' @seealso \code{\link{plot.xbal}}, \code{\link{xBalance}},
#' \code{\link{segments}}, \code{\link{points}}
#' @example inst/examples/balanceplot.R
#' @export
#' @import grDevices
balanceplot <- function(x,
                        ordered = FALSE,
                        segments = TRUE,
                        colors = "black",
                        shapes = c(15, 16, 17, 18, 0, 1, 10, 12, 13, 14),
                        segments.args = list(col = "grey"),
                        points.args = list(cex = 1),
                        xlab = "Balance",
                        xrange = NULL,
                        groups = NULL,
                        tiptext = NULL,
                        include.legend = TRUE,
                        legend.title = NULL,
                        ...) {

  stopifnot(length( dx <- dim(x) ) == 2, dx >= 1)
  names(dx) <- NULL # for comparisons below
  nvars  <- dx[1]
  nstrat <- dx[2]

  if (is.null(rownames(x))) {
    rownames(x) <- paste0("V", 1:nvars)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("S", 1:nstrat)
  }

  # create some default tooltips if needed, will only be used if user wraps this in RSVGTipsDevice
  if (is.null(tiptext)) {
    tiptext <- array(data = c(rep(rownames(x), dx[2]),
                              rep(colnames(x), each = dx[1])),
                     c(dx, 2))

  }

  # just make sure that colors and shapes have the right length/shape
  if (is.vector(colors)) {
    origcolors<-colors
    colors <- rep(colors, length.out = nstrat)
    colors <- matrix(rep(colors, each = nvars), nrow = nvars)
  }

  if (!identical(dim(colors), dx)) {
    stop("`colors` argument must have the same dims as `x`, or be comformable.")
  }

  if (is.vector(shapes)) {
    origshapes<-shapes
    shapes <- rep(shapes, length.out = nstrat)
    shapes <- matrix(rep(shapes, each = nvars), nrow = nvars)
  }

  if (!identical(dim(shapes), dx)) {
    stop("`shapes` argument must have the same dims as `x`, or be comformable.")
  }

  ngrps <- 0
  if (!is.null(groups)) {
    nagrp <- is.na(groups)
    ngrps <- length(unique(groups[which(!nagrp)]))
  }

  if(is.null(xrange)){
  xrange <- range(x, na.rm = TRUE)
  xrange <- xrange + xrange * 0.25
  ##xrange <- range(x, na.rm = TRUE) * 1.25
  }
  # we want a line for each of the variables, two lines for each group, and extra lines for the legend equal to the number of stata, and one for the optional legend title.
  yrange <- c(1, nvars + 2 * ngrps + 1 + nstrat + ifelse(!is.null(legend.title), 1, 0))

  if (ordered) {
    # order X by the groups, and within groups order by the first column
    localorder <- order( x[,1])
    x <- x[localorder, , drop = F]
  }

  if (!is.null(groups)) {
    rownames(x) <- paste0(rownames(x),
                          ifelse(is.na(groups), "", "    "))
  }

  if (names(dev.cur()) != "svg") {
    mai <- par('mai')
    mai[2] <- max(strwidth(rownames(x), units = "inches")) + mai[2]
    par(mai = mai)
  } else {
    mar <- par("mar")
    mar[2] <- max(nchar(x)) + mar[2] # assume one line per character
    par(mar = mar)
  }

  plot(xrange,
       yrange,
       axes = FALSE,
       pch = 19,
       col = "blue",
       ylab = "",
       xlab = xlab,
       type = "n",
       ...)

  axis(1, at = pretty(seq(xrange[1], xrange[2], length = 5)))

  if (is.null(groups)) {

    .balanceplot(x, segments, shapes, colors, segments.args, points.args, 0, tiptext)

  } else {
    offset <- 0
    nagrp <- is.na(groups)
    gnames <- unique(na.omit(groups))
    for (g in gnames) {

      idx <- groups == g & !is.na(groups)
      subx <- x[idx,, drop = FALSE]
      subtip <- tiptext[idx,,, drop = FALSE]
      subshape <- shapes[idx,, drop = FALSE]
      subcolor <- colors[idx,, drop = FALSE]

      offset <- .balanceplot(subx, segments, subshape, subcolor, segments.args, points.args, offset, subtip)

      axis(2, labels = g, at = offset + 0.25, las = 2, tick = FALSE)

      offset <- offset + 1
    }

    if (sum(nagrp) > 0) {

      subx <- x[nagrp,, drop = FALSE]
      subtip <- tiptext[nagrp,,, drop = FALSE]
      subshape <- shapes[nagrp,, drop = FALSE]
      subcolor <- colors[nagrp,, drop = FALSE]
      .balanceplot(subx, segments, subshape, subcolor, segments.args, points.args, offset, subtip)
    }

  }

  abline(v = 0, col = "#333333")


  if (length(colnames(x)) > 0 && include.legend) {
    legend(x = "topright",
           legend = colnames(x),
           pch = origshapes,
           col = origcolors,
           title = legend.title,
           bty = "n")
  }

}

.balanceplot <- function(x, segments, shapes, colors, segments.args, points.args, offset, tiptext) {
  n <- dim(x)[1]
  nstrat <- dim(x)[2]
  ypos <- n:1 + offset

  tts <- "devSVG" == names(dev.cur())[1] && requireNamespace("RSVGTipsDevice")

  if (segments && dim(x)[2] > 1) {
    bnds <- t(apply(x, 1, range))
    do.call(graphics::segments,
            append(list(x0 = bnds[,1],
                        y0 = ypos,
                        x1 = bnds[,2],
                        y1 = ypos),
                   segments.args))
  }

  for(i in 1:nstrat) {

    for (j in seq_along(ypos)) {

      if (tts) {
        # note that these indices are reversed versus convention [i, j, k] notation
        # i is strata (the columns of our tiptext object)
        # j is the variable (the rows of the tips)
        if (dim(tiptext)[3] == 2) {
          RSVGTipsDevice::setSVGShapeToolTip(tiptext[j, i, 1], tiptext[j, i, 2])
        }
        if (dim(tiptext)[3] == 1) {
          RSVGTipsDevice::setSVGShapeToolTip(tiptext[j, i, 1])
        }
      }

      do.call(graphics::points,
              append(list(x[j, i],
                          ypos[j],
                          pch = shapes[j, i],
                          col = colors[j, i]),
                     points.args))
    }
  }


  axis(2, labels = rownames(x), at = ypos, las = 2, tick = FALSE)

  return(offset + n + 1)
}
