.hst.main <- 
function(x, col.fill, col.stroke, col.bg, col.grid,
       col.box, col.reg,
       over.grid, cex.axis, col.axis, rotate.values, offset,
       breaks, bin.start, bin.width,
       bin.end, prop, hist.counts, cumul,
       xlab, ylab, main, sub, quiet, fun.call=NULL, ...) {


  # scale for regular R or RStudio
  adj <- .RSadj(bubble.size=NULL, cex.axis)
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab

  if (is.numeric(breaks) && !is.null(bin.start)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Choose only one option to specify a start value.\n",
      "Either choose the option  breaks  or the option  bin.start.\n\n")
  }

  # get variable labels if exist plus axes labels
  if (is.null(ylab))
    ylab <- ifelse (!prop, "Count of", "Proportion of")
  gl <- .getlabels(xlab, ylab, main, sub, cex.lab=getOption("lab.size"))
  x.name <- gl$xn; x.lbl <- gl$xl
  x.lab <- gl$xb
  y.lab <- paste(gl$yb, x.name)
  main.lab <- gl$mb
  sub.lab <- gl$sb
  cex.lab <- gl$cex.lab

  # get breaks from user supplied bin width and/or supplied start value
  if (!is.null(bin.width)  || !is.null(bin.start) || !is.null(bin.end)) {
    if (is.null(bin.start)) 
      bin.start <- pretty(min(x, na.rm = TRUE):max(x, na.rm = TRUE))[1]
    if (is.null(bin.width)) {
      h <- suppressWarnings(hist(x, plot=FALSE, breaks="Sturges"))
      bin.width <- h$breaks[2]-h$breaks[1]
    }
    max.x <- max(x, na.rm=TRUE)
    if (is.null(bin.end)) bin.end <- max.x
    if (bin.end < bin.start) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "bin.start: ", bin.start, "\n",
        "bin.end: ", bin.end, "\n",
        "bin.end is larger than bin.start, make bin.end larger.\n\n")
    }
    breaks <- seq(bin.start,bin.end,bin.width)
    seq.end <- bin.end
    while (max(breaks) < bin.end) {
      seq.end <- seq.end + bin.width
      breaks <- seq(bin.start,seq.end,bin.width)
    }
  }
  
  # for user supplied bins, from seq function or bin.start, 
  #  make sure entire data range is spanned
  if (is.numeric(breaks)) {
    cc <- cut(x, breaks, dig.lab=6, ...)   # replace each data value with its bin
    labs <- levels(cc)  # get list of unique bins, ordered
    bins <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
          upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    bin.min <- min(bins)
    bin.max <- max(bins)
    n.under <- length(x[which(x<bin.min)])
    n.over <- length(x[which(x>bin.max)])
    if (n.under+n.over > 0) {
      txt.u <- "";  txt.o <- "";  txt.nu <- "";  txt.no <- ""
      if (length(breaks) > 3)
        txt.c <- paste("Specified bin cutpoints: ", bin.min, breaks[2], "...", 
          breaks[length(breaks)-1], bin.max, "\n\n")
      else
        txt.c <- paste("Range of the specified bins: ", bin.min,
          " to ", bin.max, "\n", sep="")
      if (n.under > 0) 
        txt.u <- paste("Data values too small to fit in the bins: ",
          x[which(x<bin.min)], "\n\n")
      if (n.over > 0)
        txt.o <- paste("Data values too large to fit in the bins: ",
          x[which(x>bin.max)], "\n\n")
      txt <- "To fix this problem, extend the bin range "
      if (n.under > 0)
        txt.nu <- paste(txt, "below ", bin.min, "\n", sep="")
      if (n.over > 0)
        txt.no <- paste(txt, "above ", bin.max, "\n", sep="")
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Range of the data for ", x.name, ": ", min(x, na.rm=TRUE), " to ",
            max(x, na.rm=TRUE), "\n",
        txt.c,
        txt.u,
        txt.o,
        "Each data value must be in a bin\n",
        txt.nu,
        txt.no, "\n",
        "Extend the bin range by setting bin.start and rerun\n\n")
    }
  }


  # calculate but do not plot the histogram
  # arguments in ... for plotting instructions generate warnings with no plot
  h <- suppressWarnings(hist(x, plot=FALSE, breaks, labels=hist.counts, ...))
  
  # relative frequency histogram option
  if (prop) h$counts <- h$counts/length(x)
    
  # cumulative histogram option
  if (cumul != "off") {
    old.counts <- h$counts
    h$counts <- cumsum(h$counts)
  }

  if (is.null(main)) {
    orig.params <- par(no.readonly=TRUE)
    on.exit(par(orig.params))
    par(mar=c(4,4,2,2)+0.1)
  }
  
  # set up plot area
  plot(h, freq=TRUE, axes=FALSE, ann=FALSE, ...)

  # axis, axis ticks
  .axes(x.lvl=NULL, y.lvl=NULL, axTicks(1), axTicks(2),
        par("usr")[1], par("usr")[3], size.axis, col.axis,
        rotate.values, offset, ...)

  # axis labels
  max.lbl <- max(nchar(axTicks(2)))
  .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl, 
          xy.ticks=TRUE, offset=offset, cex.lab=size.lab, ...) 

  # colored background for plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border=col.box)
  
  # grid lines computation
  vy <- pretty(h$counts)
  vx <- h$breaks

  # plot the histogram and grid lines
  if (!over.grid) {
    abline(v=seq(vx[1],vx[length(vx)],vx[2]-vx[1]), col=col.grid, lwd=.5)
    abline(h=seq(vy[1],vy[length(vy)],vy[2]-vy[1]), col=col.grid, lwd=.5)
  }

  plot(h, add=TRUE, col=col.fill, border=col.stroke, freq=TRUE, labels=hist.counts, ...)
  if (cumul == "both") {
    h$counts <- old.counts
    plot(h, add=TRUE, col=col.reg, freq=TRUE)
  }
  if (over.grid) {
    abline(v=seq(vx[1],vx[length(vx)],vx[2]-vx[1]), col=col.grid, lwd=.5)
    abline(h=seq(vy[1],vy[length(vy)],vy[2]-vy[1]), col=col.grid, lwd=.5)
  }

 
#------------
# text output
#------------
  if (!quiet) {

    stats <- .hst.stats(h, length(x), fun.call)

    txsug=stats$txsug
    tx=stats$tx
    bin.width=stats$bin.width
    n.bins=stats$n.bins
    prop=stats$prop
    cum.c=stats$counts_cum
    cum.p=stats$prop_cum

    return(list(txsug=txsug, ttx=tx, bin.width=bin.width, n.bins=n.bins, breaks=h$breaks, 
      mids=h$mids, counts=h$counts, prop=prop, counts_cum=cum.c, prop_cum=cum.p))
  }

}
