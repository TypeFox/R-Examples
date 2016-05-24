.dn.main <- 
function(x, 
         bw, type, bin.start, bin.width,
         col.fill, col.bg, col.grid, col.box, col.nrm, col.gen,
         col.fill.nrm, col.fill.gen,
         cex.axis, col.axis, rotate.values, offset, 
         x.pt, xlab, main, sub, y.axis, x.min, x.max, band, quiet, ...)  {


  # scale for regular R or RStudio
  adj <- .RSadj(bubble.size=NULL, cex.axis)
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab

  if (!is.null(x.pt)) {
    y.axis <- TRUE
    type <- "general"
  }

  # get variable labels if exist plus axes labels
  gl <- .getlabels(xlab, main=main, cex.lab=getOption("lab.size"))
  x.name <- gl$xn; x.lbl <- gl$xl;
  x.lab <- gl$xb
  main.lab <- gl$mb
  sub.lab <- gl$sb
  cex.lab <- gl$cex.lab

  # get breaks from user supplied bin width and/or supplied start value
  # otherwise, breaks="Sturges by default
  # copied from histogram
  if (!is.null(bin.width)  || !is.null(bin.start)) {
    if (is.null(bin.start)) bin.start <- pretty(min(x):max(x))[1]
    if (is.null(bin.width)) {
      h <- hist(x, plot=FALSE, breaks="Sturges")
      bin.width <- h$breaks[2]-h$breaks[1]
    }
    max.x <- max(x, na.rm = TRUE)
    seq.end <- max.x
    breaks <- seq(bin.start,seq.end,bin.width)
    while (max(breaks) < max.x) {
      seq.end <- seq.end + bin.width
      breaks <- seq(bin.start,seq.end,bin.width)
    }
  }
  else breaks="Sturges"


  if (is.null(main)) {
    orig.params <- par(no.readonly=TRUE)
    on.exit(par(orig.params))
    par(mar=c(4,4,2,2)+0.1)
  }

  # histogram calculations, no plot
  h <- hist(x, plot=FALSE, breaks)

  n <- sum(!is.na(x))
  n.miss <- sum(is.na(x))
  if (n.miss > 0) x <- na.omit(x)
 
  # general density curve, no plot
  # suppress warnings about possible graphic parameters
  d.gen <- suppressWarnings(density(x, bw, ...))
  
  mx <- mean(x)

  # min and max x coordinates for graph, make symmetric
  min.dev.x <- min(d.gen$x) - mx
  max.dev.x <- max(d.gen$x) - mx
  if (abs(min.dev.x) > abs(max.dev.x)) {
    if (is.null(x.min)) x.min <- min(d.gen$x)
    if (is.null(x.max)) x.max <- mx + abs(min.dev.x)
  }
  if (abs(max.dev.x) >= abs(min.dev.x)) {
    if (is.null(x.min)) x.min <- mx - abs(max.dev.x)
    if (is.null(x.max)) x.max <- max(d.gen$x)
  }
  
  # normal density curve, no plot
  xx <- seq(x.min, x.max, length=200)
  if (col.fill.nrm == "transparent") lw <- 2 else lw <- 1
  d.nrm <- dnorm(xx,mean(x),sd(x))

  # max y coordinate for graph
  max.y <- max(max(d.nrm), max(d.gen$y), max(h$density))

  # set up plot area
  # bw if specified also gets passed to plot, so suppress warning

  plot(h, border="transparent", freq=FALSE,
     xlim=c(x.min,x.max), ylim=c(0,max.y),
     axes=FALSE, ann=FALSE, xlab=NULL, ylab=NULL, main=NULL, ...)

  # axis, axis ticks
  if (!y.axis)
    .axes(x.lvl=NULL, y.lvl=NULL, axTicks(1), NULL,
          par("usr")[1], par("usr")[3], size.axis, col.axis,
          rotate.values, offset, ...)
  else
    .axes(x.lvl=NULL, y.lvl=NULL, axTicks(1), axTicks(2),
          par("usr")[1], par("usr")[3], size.axis, col.axis,
          rotate.values, offset, ...)

  # axis value labels
  if (!y.axis) y.lab="" else y.lab="Density"
  max.lbl <- max(nchar(axTicks(2)))
  .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl, 
          xy.ticks=TRUE, offset=offset, cex.lab=size.lab, ...) 

  # colored background for plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border=col.box)
  
  # plot the histogram
  plot(h, add=TRUE, freq=FALSE, col=col.fill, border="transparent")

  # plot the normal curve
  if (type == "normal" || type == "both") {
    if (col.fill.nrm == "transparent") lw <- 1.35 else lw <- 1
    polygon(c(x.min,xx,x.max), c(0,d.nrm,0), col=col.fill.nrm, 
            border=col.nrm, lwd=lw)
  }

  # plot the general curve
  if (type == "general" || type == "both") {
    if (col.fill.gen == "transparent") lw <- 1.35 else lw <- 1
    polygon(d.gen, col=col.fill.gen, border=col.gen, lwd=lw)
  }

  if (band) rug(x, col="gray70")
 
  # plot the optional bar about a chosen point for general curve only
  if (!is.null(x.pt)  &&  type == "general") {
    d <- d.gen
    y.pt <- d$y[which.min(abs(x.pt-d$x))]
    xbeg <- x.pt - 0.5
    xend <- x.pt + 0.5
    xsub <- d$x[d$x > xbeg & d$x < xend]
    ysub <- d$y[d$x > xbeg & d$x < xend]
    txt <- paste("Density =", .fmt(y.pt,3), "at x =", .fmt(x.pt, 2))
    title(main=txt)
    polygon(c(xbeg, xsub, xend), c(0, ysub, 0), col="lightsteelblue")
    if (min(x) > 0) left <- -2*min(x) else left <- 2*min(x) 
    lines(c(left,x.pt), c(y.pt,y.pt), col="darkblue", lwd = 0.5)
    lines(c(x.pt,x.pt), c(0,y.pt), col="darkblue", lwd = 0.5)
  }

  # text output
  tx=""
  if (!quiet) {

    tx <- character(length = 0)

    tx[length(tx)+1] <- paste("Sample Size: ", n)
    tx[length(tx)+1] <- paste("Missing Values: ", n.miss)
 
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("Density bandwidth for general curve: ",
      .fmt(d.gen$bw,4), sep="")
    tx[length(tx)+1] <- "For a smoother curve, increase bandwidth with option: bw"

    tx[length(tx)+1] <- ""
    W <- NA; p.val <- NA
    if (type == "normal" || type == "both") {
      digits.d <- 3
      if (n > 2 && n < 5000) {
        nrm <- shapiro.test(x)
        W <- .fmt(nrm$statistic,min(4,digits.d+1))
        p.val <- .fmt(nrm$p.value,min(4,digits.d+1))
        tx[length(tx)+1] <- paste("Null hypothesis is a normal population")
        tx[length(tx)+1] <- paste(nrm$method, ":  W = ", W, ",  p-value = ",
          p.val, sep="")
      }
      else 
        tx[length(tx)+1] <- "Sample size out of range for Shapiro-Wilk normality test."
    }

  return(list(tx=tx, bw=d.gen$bw, n=n, n.miss=n.miss, W=W, pvalue=p.val))

  }

} 
