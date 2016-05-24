## likert <- function(x, ..., xName=deparse(substitute(x))) ## simplifies calling sequence
##   plot.likert(x, ..., xName=xName)

if.R(r={
  plot.likert <- function(x, ...)
    UseMethod("plot.likert")
}, s={
  plot.likert <- function(x, ...)
    stop("The likert functions are not available in S-Plus.")
})

likert <- plot.likert
likertplot <- plot.likert

xscale.components.top.HH <- function(...) {
  ans <- xscale.components.default(...)
  ans$top <- ans$bottom
  ans$top$labels$labels <- names(ans$top$labels$labels) ## this requires named bottom labels!
  ans
}
## environment(xscale.components.top.HH) <- environment(plot.likert)

yscale.components.right.HH <- function(...) {
  ans <- yscale.components.default(...)
  ans$right <- ans$left
  ans$right$labels$labels <- names(ans$right$labels$labels) ## this requires named left labels!
  ans
}
## environment(yscale.components.right.HH) <- environment(plot.likert)


## panel.likert <- function(...)
##   .Defunct("panel.barchart", package="HH")
## ## panel.likert <- function(..., rightAxisLabels, rightAxis) {
## ##   panel.barchart(...)
## ##   if (rightAxis) panel.axis.right(side="right", at=1:length(rightAxisLabels),
## ##                                   labels=rightAxisLabels, outside=TRUE)
## ## }

plot.likert.default <- function(x,
                                positive.order=FALSE,
                                ylab=names(dimnames(x)[1]),
                                xlab=if (as.percent != FALSE) "Percent" else "Count",
                                main=xName,
                                reference.line.col="gray65",
                                col.strip.background="gray97",
                                col=likertColor(attr(x, "nlevels"),
                                  ReferenceZero=ReferenceZero,
                                  colorFunction=colorFunction,
                                  colorFunctionOption=colorFunctionOption),
                                colorFunction="diverge_hcl",
                                colorFunctionOption="lighter",
                                as.percent=FALSE,
                                par.settings.in=NULL,
                                horizontal=TRUE,
                                ReferenceZero=NULL,
                                ...,
                                key.border.white=TRUE,
                                xName=deparse(substitute(x)),
                                rightAxisLabels=rowSums(abs(x)),
                                rightAxis=!missing(rightAxisLabels),
                                ylab.right=if (rightAxis) "Row Count Totals" else NULL,
                                panel=panel.barchart,
                                xscale.components=xscale.components.top.HH,
                                yscale.components=yscale.components.right.HH,
                                xlimEqualLeftRight=FALSE,
                                xTickLabelsPositive=TRUE,
                                reverse=FALSE) {
  force(xName)
  rightAxisMissing <- missing(rightAxis)  ## needed by as.percent
  x.input <- x
  if (is.null(dim(x))) {
    x <- t(x)
    if (is.null(dimnames(x))) dimnames(x) <- list("", letters[seq(along=x)])
    dimnames(x)[[1]] <- ""
  }
  force(rightAxis)
  force(rightAxisLabels)
  force(ylab.right)
  if (as.percent != FALSE) {
    x.pct <- x / rowSums(abs(x)) * 100
    x.pct[x==0] <- 0
    x <- as.likert(x.pct,
                   ReferenceZero=ReferenceZero)
    if (rightAxisMissing && as.percent != "noRightAxis" ) {
      rightAxis <- TRUE
      if (is.null(ylab.right))
        ylab.right <- "Row Count Totals"
    }
    ## else
    ##   rightAxis <- FALSE
  } else {
    x <- as.likert(x,
                   ReferenceZero=ReferenceZero)
  }

  if (!is.null(ReferenceZero) && !is.null(attr(x, "ReferenceZero"))) {
    if (ReferenceZero != attr(x, "ReferenceZero"))
      warning(paste('(Argument ReferenceZero = ', ReferenceZero, ') != (',
                    'as.likert ReferenceZero = ', attr(x, "ReferenceZero"), ')\n',
                    'as.likert ReferenceZero =', attr(x, "ReferenceZero"),
                    'will be used.'))
  }
  if (is.null(ReferenceZero) && !is.null(attr(x, "ReferenceZero")))
    ReferenceZero <- attr(x, "ReferenceZero")
  if (!is.null(ReferenceZero) && is.null(attr(x, "ReferenceZero")))
    warning(paste('(Argument ReferenceZero = ', ReferenceZero, ') != (',
                  'as.likert ReferenceZero = NULL)\n',
                  'Argument ReferenceZero will be ignored.'))

  key.title <- names(dimnames(x)[2])
  if (is.null(key.title)) key.title <- attr(x,"names.dimnames")[[2]]
  if (is.null(key.title)) key.title <- " "
  auto.key.likert <- list(title=key.title,
                          lines.title=1.5,
                          text=attr(x, "original.levels"),
                          cex=.7,
                          border=FALSE,
                          height=1,
                          space="bottom",
                          columns=attr(x, "nlevels"),
##                          columns=min(2, length(attr(x, "original.levels"))), ## attr(x, "nlevels"),
                          padding.text=1,
                          size=2,
                          between=.5,
                          between.columns=2,
                          just=.5,
                          reverse=FALSE,
                          rect=list(col=col, border=if (key.border.white) "white" else col),
                          ##                ## The next two lines suppress unwanted automatic displays.
                          points=FALSE,     ## This line is necessary when the right axis is used.
                          rectangles=FALSE) ## This line is necessary and not redundant.

  dotdotdot <- list(...)
  if (!is.null(dotdotdot$auto.key)) {
    ak <- dotdotdot$auto.key
    auto.key.likert[names(ak)] <- ak
    dotdotdot$auto.key <- NULL
  }
  ## auto.key.likert$rect=list(
  ##   col=col, border=col,
  ##   height=auto.key.likert$height, size=auto.key.likert$size)

  dotdotdot$scales$alternating <- 1

  if (missing(ylab) && (is.null(ylab)||is.na(ylab))) ylab <- NULL
  ## RColorBrewer diverging palettes: c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  ## These are the middle colors from RCOlorBrewer:
  ## > for (i in c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdGy", "RdYlBu", "RdYlGn", "Spectral"))
  ## + print(c(i, RColorBrewer::brewer.pal(n=3, name=i)[2]))
  ## [1] "RdBu"     "#F7F7F7"
  ## [1] "BrBG"     "#F5F5F5"
  ## [1] "PiYG"     "#F7F7F7"
  ## [1] "PRGn"     "#F7F7F7"
  ## [1] "PuOr"     "#F7F7F7"
  ## [1] "RdGy"     "#FFFFFF"
  ## [1] "RdYlBu"   "#FFFFBF"
  ## [1] "RdYlGn"   "#FFFFBF"
  ## [1] "Spectral" "#FFFFBF"

  ## nc <- ncol(x)
  ## if (missing(middle.color)) middle.color ## "#F7F7F7"
  ## ## if middle.color is missing as an argument, then use the default value from the argument list
  ## if (attr(x, "even.col")) { ## no zero to split
  ##   likert.palette <- col[c((nc/2):1, ((nc/2)+1):nc)]
  ## }
  ## else { ## yes zero to split
  ##   likert.palette <- col[c((nc/2):1, ((nc/2)+1):(nc-1))]
  ## }
  likert.palette <- col[attr(x, "color.seq")]

  if (positive.order) {
    ## x.attr <- attributes(x)
    ## x.attr$dim <- NULL
    ## x.attr$dimnames <- NULL
    ## x.attr$original.order <- order(x.attr$positive.order)
    ## x[] <- x[x.attr$positive.order,, drop=FALSE]
    ## attributes(x)[names(x.attr)] <- x.attr
    ## rightAxisLabels <- rightAxisLabels[x.attr$positive.order] ## rev(rev(rightAxisLabels)[x.attr$positive.order])
    pos.order <- attr(x, "positive.order")
    attr(x, "original.order") <- order(pos.order)
    x[] <- x[pos.order, , drop=FALSE]
    dimnames(x)[[1]] <- dimnames(x)[[1]][pos.order]
    rightAxisLabels <- rightAxisLabels[pos.order]
  }
  if ((horizontal + positive.order + reverse) %% 2) { ## if one or three, then reverse
    x <- rev(x)
    rightAxisLabels <- rev(rightAxisLabels)
  }

  par.settings <- list(strip.background=list(col=col.strip.background),
                       reference.line=list(col=reference.line.col),
                       layout.heights=list(
                         main.key.padding=2.5,
                         key.axis.padding=0,
                         axis.top=.75,
                         xlab.key.padding=2),
                       layout.widths=list(
                         ylab.right=if (rightAxis) 5 else
                         trellis.par.get("layout.widths")$ylab.right,
                         right.padding=if (rightAxis) 0 else
                         trellis.par.get("layout.widths")$right.padding),
                       clip=list(panel="off"))
  par.settings[names(par.settings.in)] <- par.settings.in
  barchart.args <- list(x=x,
                        as.table=TRUE,
                        col=likert.palette,
                        border=likert.palette,
                        auto.key=auto.key.likert,
                        xlab=xlab, ylab=ylab,
                        ylab.right=ylab.right,
                        par.settings=par.settings,
                        reference.line=TRUE,
                        main=main,
                        horizontal=horizontal,
                        ## rightAxisLabels=rightAxisLabels,
                        ## rightAxis=rightAxis,
                        panel=panel,
                        xscale.components=xscale.components,
                        yscale.components=yscale.components)
  barchart.args[names(dotdotdot)] <- dotdotdot
  ## if (!is.null(barchart.args$horizontal) && !barchart.args$horizontal) {
  ##   tmp <- barchart.args$xlab
  ##   barchart.args$xlab <- barchart.args$ylab
  ##   barchart.args$ylab <- tmp
  ## } else {
  ##   barchart.args$horizontal <- TRUE
  ## }
  if (is.null(barchart.args$horizontal))
    barchart.args$horizontal <- TRUE
  if (!barchart.args$horizontal) {
    tmp <- barchart.args$xlab
    barchart.args$xlab <- barchart.args$ylab
    barchart.args$ylab <- tmp
  }
  result <- do.call("barchart", barchart.args)
  if (rightAxis) {
    if (barchart.args$horizontal) {
      result$y.scales$alternating <- 3
      names(result$y.limits) <- rightAxisLabels ## rev(rightAxisLabels)
      result$y.scales$tck <- c(1, 1)
      result$y.scales$col.line <- 0
    } else {
      result$x.scales$alternating <- 3
      names(result$x.limits) <- rightAxisLabels
    }
    ## class(result) <- c("trellis.right.HH", class(result))
  }

  if (horizontal) {
    if (xlimEqualLeftRight
        ## &&
        ## (is.null(list(...)$xlim) &&
        ##  is.null(list(...)$scales$limits) &&
        ##  is.null(list(...)$scales$x$limits))
        ) result$x.limits <- c(-1,1) * max(abs(result$x.limits))

    if (xTickLabelsPositive
        &&
        (is.null(list(...)$scales$labels) &&
         is.null(list(...)$scales$x$labels))
        ) {
      if (!is.numeric(result$x.scales$at)) {
        x.range <- result$x.limits
        tn <- list(...)$scales$x$tick.number
        if (is.null(tn)) tn <- list(...)$scales$tick.number
        if (is.null(tn)) tn <- formals(pretty.default)$n
        result$x.scales$at <- pretty(x.range, tn)
      }
      result$x.scales$labels <- abs(result$x.scales$at)
    }
  } else { ## vertical
    if (xlimEqualLeftRight
        ## &&
        ## (is.null(list(...)$ylim) &&
        ##  is.null(list(...)$scales$limits) &&
        ##  is.null(list(...)$scales$y$limits))
        ) result$y.limits <- c(-1,1) * max(abs(result$y.limits))

    if (xTickLabelsPositive
        &&
        (is.null(list(...)$scales$labels) &&
         is.null(list(...)$scales$y$labels))
        ) {
      if (!is.numeric(result$y.scales$at)) {
        y.range <- result$y.limits
        tn <- list(...)$scales$y$tick.number
        if (is.null(tn)) tn <- list(...)$scales$tick.number
        if (is.null(tn)) tn <- formals(pretty.default)$n
        result$y.scales$at <- pretty(y.range)
      }
      result$y.scales$labels <- abs(result$y.scales$at)
    }
  }

  result$axis <- axis.RightAdjustRight
  result
}
## environment(plot.likert.default) <- environment(plot.likert)
## assignInNamespace("plot.likert.default", plot.likert.default, "HH")
plot.likert.array <- function(x,  ## an array
                              condlevelsName=paste("names(dimnames(", xName, "))[-(1:2)]", sep=""),
                              xName=deparse(substitute(x)),
                              main=paste("layers of", xName, "by", condlevelsName),
                              ...) {
  ##force(condlevelsName)
  force(xName)
  if (length(dim(x))==2) NextMethod("plot.likert")
  x <- as.MatrixList(x)  ## list of matrices, one per each layer of array
  plot.likert.list(x,
                   condlevelsName=condlevelsName,
                   xName=xName,
                   main=main,
                   ...)
}
  ## x <- lapply(x, as.likert, ...)
  ## x.pl <- lapply(x, plot.likert, positive.order=positive.order, horizontal=horizontal, ...) ## named list of likert plots
  ## x.pl.nonames <- x.pl ## if (strip.left) about to become unnamed list of likert plots
  ## names(x.pl.nonames) <- NULL ## names are removed

  ## if (strip.left) {
  ##   ResizeEtc.likertPlot(do.call("c", x.pl.nonames),
  ##                        x,
  ##                        x.pl.nonames,
  ##                        horizontal=horizontal,
  ##                        condlevelsName=condlevelsName,
  ##                        x.same=horizontal, y.same=!horizontal,
  ##                        layout=layout,
  ##                        strip=strip,
  ##                        strip.left.values=strip.left.values,
  ##                        strip.left.par=strip.left.par,
  ##                        resize.height=resize.height,
  ##                        resize.width=resize.width,
  ##                        main=main)
  ## } else {
  ##   ResizeEtc.likertPlot(do.call("c", x.pl.nonames),
  ##                        x,
  ##                        x.pl.nonames,
  ##                        horizontal=horizontal,
  ##                        condlevelsName=condlevelsName,
  ##                        x.same=horizontal, y.same=!horizontal,
  ##                        layout=layout,
  ##                        strip=strip,
  ##                        strip.values=strip.values,
  ##                        resize.height=resize.height,
  ##                        resize.width=resize.width,
  ##                        main=main)
  ## }
  ## }
if (FALSE) {
plot.likert.array <- function(x,  ## an array
                              condlevelsName=paste(names(dimnames(x))[-(1:2)], collapse="."),
                              xName=deparse(substitute(x)),
                              main=paste("layers of", xName, "by", condlevelsName),
                              layout=if (horizontal) c(1, length(dim(x))-2) else c(length(dim(x))-2, 1),
                              positive.order=FALSE,
                              strip=TRUE,
                              strip.left=TRUE,
                              strip.values=names(x.pl),
                              strip.left.par=list(cex=1, lines=1),
                              horizontal=TRUE,
                              ...,
                              resize.height=c("nrow","rowSums"),
                              resize.width=1) {
  force(condlevelsName)
  force(xName)
  if (length(dim(x))==2) NextMethod("plot.likert")
  if (class(resize.height)=="character") {
    resize.height <- switch(match.arg(resize.height),
                            nrow=rep(dim(x)[1], layout[2])+1,
                            rowSums=apply(x, length(dim(x)), function(x) sum(abs(x))),
                            stop("invalid value for resize.height"))
  }
  if (length(resize.height) != layout[2])
    stop("Wrong length for resize.height.")
  x <- as.MatrixList(x)  ## list of matrices, one per each layer of array
  x <- lapply(x, as.likert, ...)
  x.pl <- lapply(x, plot.likert, positive.order=positive.order, horizontal=horizontal, ...) ## named list of likert plots
  x.pl.nonames <- x.pl ## if (strip.left) about to become unnamed list of likert plots
  names(x.pl.nonames) <- NULL ## names are removed

  if (strip.left) {
    ResizeEtc.likertPlot(do.call("c", x.pl.nonames),
                         x,
                         x.pl.nonames,
                         horizontal=horizontal,
                         condlevelsName=condlevelsName,
                         x.same=horizontal, y.same=!horizontal,
                         layout=layout,
                         strip=strip,
                         strip.left.values=strip.left.values,
                         strip.left.par=strip.left.par,
                         resize.height=resize.height,
                         resize.width=resize.width,
                         main=main)
  } else {
    ResizeEtc.likertPlot(do.call("c", x.pl.nonames),
                         x,
                         x.pl.nonames,
                         horizontal=horizontal,
                         condlevelsName=condlevelsName,
                         x.same=horizontal, y.same=!horizontal,
                         layout=layout,
                         strip=strip,
                         strip.values=strip.values,
                         resize.height=resize.height,
                         resize.width=resize.width,
                         main=main)
  }
}
## environment(plot.likert.array) <- environment(plot.likert)
}

plot.likert.list <- function(x,  ## named list of matrices, 2D tables, 2D ftables, or 2D structables, or all-numeric data.frames
                             condlevelsName="ListNames",
                             xName=deparse(substitute(x)),
                             main=paste("List items of", xName, "by", condlevelsName),
                             layout=if (length(dim.x) > 1) dim.x else {
                               if (horizontal) c(1, length(x)) else c(length(x), 1)},
                             positive.order=FALSE,
                             strip=!horizontal,
                             strip.left=horizontal,
                             strip.left.values=names(x),
                             strip.values=names(x),
                             strip.par=list(cex=1, lines=1),
                             strip.left.par=list(cex=1, lines=1),
                             horizontal=TRUE,
                             ...,
                             rightAxisLabels=sapply(x, function(x) rowSums(abs(x)), simplify=FALSE),
                             rightAxis=!missing(rightAxisLabels),
                             resize.height.tuning=-.5,
                             resize.height=if (missing(layout) || length(dim.x) != 2) {
                               c("nrow","rowSums")
                             } else {
                               rep(1, layout[2])
                             },
                             resize.width=if (missing(layout)) {
                               1
                             } else {
                               rep(1, layout[1])
                             },
                             box.ratio=if (
                               length(resize.height)==1 &&
                               resize.height == "rowSums") 1000 else 2,
                             xscale.components=xscale.components.top.HH,
                             yscale.components=yscale.components.right.HH) {
  force(xName)
  ##force(layout)
  ##force(resize.height)
  ##force(resize.width)
  ##force(box.ratio)

  ## if (!is.null(dim(x))) stop(paste(xName, " has dimension=", deparse(dim(x)),
  ##                                  ". plot.likert.list requires a list without a dim attribute.", sep=""))
  for (nxi in names(x)) { ## convert vectors to single-row matrices
    xi <- x[[nxi]]
    if (is.numeric(xi) && is.null(dim(xi))) x[[nxi]] <- t(xi)
  }
  if (!is.listOfNamedMatrices(x)) {
    if (is.null(names(x)) || any(is.na(names(x)))) stop("Items in a list for plot.likert must be named.")
    if (!all(sapply(x, function(x) length(dim(x))) == 2))
      stop("All items in a list for plot.likert must have at most two dimensions.")
    if (!all(sapply(x, ncol) == ncol(x[[1]])))
      stop("All items in a list for plot.likert must have the same number of columns.")
    if (is.data.frame(x))
      stop("plot.likert.list does not accept a data.frame.\nPlease use plot.likert.data.frame.")
    ## if (any(sapply(x, function(xx) is.data.frame(xx) && !all(sapply(xx, is.numeric)))))
    ##   stop("A data.frame in a plot.likert.list argument must have only numeric columns.")
  }

  names.x <- names(x)
  dim.x <- dim(x)
  dimnames.x <- dimnames(x)
  x <- lapply(x, function(z)
              if (is.data.frame(z)) data.matrix(z[, sapply(z, is.numeric), drop=FALSE]) else z
              )
  dim(x) <- dim.x
  dimnames(x) <- dimnames.x
  names(x) <- names.x

  x.pl <- mapply(plot.likert, x,
                 rightAxisLabels=rightAxisLabels,
                 MoreArgs=list(
                   positive.order=positive.order, horizontal=horizontal, ...,
                   box.ratio=box.ratio,
                   rightAxis=rightAxis,
                   xscale.components=xscale.components,
                   yscale.components=yscale.components),
                 SIMPLIFY=FALSE, USE.NAMES=TRUE)  ## named list of likert plots
##  x.pl.nonames <- x.pl ## if (strip.left) about to become unnamed list of likert plots
##  names(x.pl.nonames) <- NULL ## names are removed

  if (class(resize.height)=="character") {
    if (resize.height=="rowSums" && !all(sapply(x, nrow)==1))
      stop("resize.height='rowSums' is not valid for panels with more than one row.")
    resize.height <- switch(match.arg(resize.height, c("nrow","rowSums")),
                            nrow=sapply(x, nrow)+resize.height.tuning,
                            rowSums=sapply(x, function(x) rowSums(abs(x)), simplify=TRUE),
                            stop("invalid value for resize.height"))
  }
  ## if (length(resize.height) != length(x))
  ##   stop("Wrong length for resize.height.")
  if (!horizontal) {
    tmp <- resize.height
    resize.height <- resize.width
    resize.width <- tmp
  }

  if (length(resize.height) > 1 && all(resize.height==resize.height[1])) resize.height <- 1
  if (length(resize.width)  > 1 && all( resize.width==resize.width[1] )) resize.width  <- 1

  if (!(length(resize.width) == 1 && length(resize.height) == 1))
    if (any(layout != c(length(resize.width), length(resize.height))))
      warning(paste("Inconsistent layout=", deparse(layout),
                    "and length(resize.width)=", deparse(length(resize.width)),
                    "and length(resize.height)=", deparse(length(resize.height))))



  result <-
    if (strip.left) {
      ResizeEtc.likertPlot(do.call("c", x.pl),
                       x,
                       x.pl,
                       horizontal=horizontal,
                       condlevelsName=condlevelsName,
                       x.same=horizontal, y.same=!horizontal,
                       layout=layout,
                       strip=strip,
                       strip.left.values=strip.left.values,
                       strip.left.par=strip.left.par,
                       resize.height=resize.height,
                       resize.width=resize.width,
                       main=main)
    } else {
  ResizeEtc.likertPlot(do.call("c", x.pl),
                       x,
                       x.pl,
                       horizontal=horizontal,
                       condlevelsName=condlevelsName,
                       x.same=horizontal, y.same=!horizontal,
                       layout=layout,
                       strip=strip,
                       strip.values=strip.values,
                       strip.par=strip.par,
                       resize.height=resize.height,
                       resize.width=resize.width,
                       main=main)
}

  if (length(dim(x)) == 2) {
    result$index.cond <- lapply(dim(x), function(i) 1:i)
    result$perm.cond <- 1:length(dim(x))
    result$condlevels <- dimnames(x)
    result$x.scales$at <- pretty(range(result$x.limits))
    result$x.scales$labels <- abs(pretty(range(result$x.limits)))
    result <- useOuterStrips(result)
  }

  result
}
## environment(plot.likert.list) <- environment(plot.likert)


ResizeEtc.likertPlot <- function(c.list,
                                 x,
                                 x.pl.nonames,
                                 horizontal,
                                 ...) {
  result <- ResizeEtc(c.list, ...)

  ## fix up axes
  if (any(unlist(lapply(x, attr, "xlimEqualLeftRight")))) {
    if (horizontal)
      result$x.limits <- c(-1, 1)*max(abs(result$x.limits))
    else
      result$y.limits <- c(-1, 1)*max(abs(result$y.limits))
  }
  if (any(unlist(lapply(x, attr, "xTickLabelsPositive")))) {
    if (horizontal) {
      xscales <- sapply(x.pl.nonames, function(x) x$x.scales[c("at","labels","tick.number")])
      winner <- which.max(sapply(xscales["at",], function(x) diff(range(x))))
      result$x.scales$at <- x.pl.nonames[[winner]]$x.scales$at
      result$x.scales$labels <- x.pl.nonames[[winner]]$x.scales$labels
    }
    else {
      yscales <- sapply(x.pl.nonames, function(x) x$y.scales[c("at","labels","tick.number")])
      winner <- which.max(sapply(yscales["at",], function(x) diff(range(x))))
      result$y.scales$at <- x.pl.nonames[[winner]]$y.scales$at
      result$y.scales$labels <- x.pl.nonames[[winner]]$y.scales$labels
    }
  }
  result
}
## environment(ResizeEtc.likert) <- environment(plot.likert)

plot.likert.table <- function(x, ..., xName=deparse(substitute(x))){
  force(xName)
  class(x) <- "array"
  plot.likert(x, xName=xName, ...)
}
plot.likert.ftable <- function(x, ..., xName=deparse(substitute(x))){
  force(xName)
  plot.likert(as.table(x), xName=xName, ...)
}
plot.likert.structable <- function(x, ..., xName=deparse(substitute(x))){
  force(xName)
  plot.likert(as.table(x), xName=xName, ...)
}
## plot.likert.numeric <- function(x, ..., xName=deparse(substitute(x))){
##   force(xName)
##   plot.likert(as.likert(x, xName=xName), xName=xName, ...)
## }
plot.likert.data.frame <- function(x, ..., xName=deparse(substitute(x))){
  force(xName)
    x.num <- data.matrix(x[, sapply(x, is.numeric), drop=FALSE]) ## not redundant, data.matrix converts character columns to NA, and factor columns to integers
  plot.likert(x.num, xName=xName, ...)
}

## The HH plot method plot.likert.likert detects "likert" objects
## created by the independent likert package and plots them correctly.
## It is not recommended that the HH package and the likert package
## both be loaded at the same time, as they have incompatible usage of
## the exported function names "likert" and "plot.likert".
plot.likert.likert <- function(x, ...) {
  ## "likert" object from independent likert package
  if (length(class(x)) == 1 && is.list(x) && !is.null(x$result))
    {
      if (is.null(x$results$Group))
        likert(Item ~ .        , data=x$results, xlab="Percent", data.order=TRUE, ...)
      else
        likert(Item ~ . | Group, data=x$results, xlab="Percent", data.order=TRUE, ...)
    }
  else
    NextMethod("plot.likert")
}



## source("c:/HOME/rmh/HH-R.package/HH/R/likert.R")
