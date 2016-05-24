plot.vgc <- function (x, y, ...,
                      m=NA, add.m=NULL, N0=NULL,
                      conf.level=.95, conf.style=c("ticks", "lines"),
                      log=c("", "x", "y", "xy"),
                      bw=zipfR.par("bw"), 
                      xlim=NULL, ylim=NULL,
                      xlab="N", ylab="V(N)", legend=NULL,
                      main="Vocabulary Growth", 
                      lty=NULL, lwd=NULL, col=NULL)
{
  ## collect all specified VGCs in single list
  VGCs <- list(x)   # this is a bit complicated because of the plot() prototype
  if (! missing(y)) {
    VGCs <- c(VGCs, list(y), list(...))
  }
  n.vgc <- length(VGCs)
  
  ## check other arguments
  log <- match.arg(log)
  conf.style <- match.arg(conf.style)
  if (!missing(m) && !missing(add.m)) stop("'m' and 'add.m' must not be specified at the same time")
  if (!missing(m) && !( is.numeric(m) && all(1 <= m & m <= 9) && length(m) == 1 ))
    stop("'m' must be a single integer between 1 and 9")
  if (!missing(add.m) && !( is.numeric(add.m) && all(1 <= add.m & add.m <= 9) ))
    stop("'add.m' must be a vector of integers between 1 and 9") 
  if (!missing(legend) && length(legend) != n.vgc)
    stop("'legend' argument must be character or expression vector of same length as number of VGCs")
  x.log <- log == "x" || log == "xy"
  y.log <- log == "y" || log == "xy"

  ## check availability of Vm in 'vgc' objects (with m or add.m option)
  req.m <- 0                            # required spectrum elements (up to req.m)
  if (!missing(m)) req.m <- m
  if (!missing(add.m)) req.m <- max(add.m)
  if (req.m > 0 && any( sapply(VGCs, function (.VGC) attr(.VGC, "m.max") < req.m) ))
    stop("all VGC arguments must include spectrum elements up to m=", req.m)
    
  ## determine maximum value that has to be accommodated on y-axis
  V.max <- max(sapply(VGCs, function (.VGC) max(V(.VGC)) ))
  if (!missing(m)) V.max <- max(sapply(VGCs, function (.VGC) max(Vm(.VGC, m)) ))
  if (!missing(add.m)) {
    for (.m in add.m) V.max <- max(V.max, sapply(VGCs, function (.VGC) max(Vm(.VGC, .m)) ))
  }
  N.max <- max(sapply(VGCs, function (.VGC) max(N(.VGC)) ))
  
  ## get default styles unless manually overridden
  if (missing(lty)) lty <- zipfR.par("lty", bw.mode=bw)
  if (missing(lwd)) lwd <- zipfR.par("lwd", bw.mode=bw)
  if (missing(col)) col <- zipfR.par("col", bw.mode=bw)

  ## typeset default label on y-axis, depending on which curves are shown
  expected <- sapply(VGCs, function (.VGC) attr(.VGC, "expected"))
  if (missing(ylab)) {
    if (!missing(m)) {                  # 'm' option => V_m(N) / E[V_m(N)]
      lab.V <- bquote(V[.(m)](N))
    }
    else {                              # default => V(N) / E[V(N)]
      lab.V <- bquote(V(N))
    }
    lab.EV <- bquote(E * group("[", .(lab.V), "]"))

    if (all(expected)) {
      ylab <- lab.EV
    }
    else if (all(!expected)) {
      ylab <- lab.V
    }
    else {
      ylab <- bquote(.(lab.V) / .(lab.EV))
    }

    if (!missing(add.m)) {              # 'add.m' option => V(N) / V_1(N) / V_2(N) / ...
      for (.m in add.m) {
        lab.Vm <- bquote(V[.(.m)](N))
        lab.EVm <- bquote(E * group("[", .(lab.Vm), "]"))
        if (any(!expected)) {
          ylab <- bquote(.(ylab) / .(lab.Vm))
        }
        if (any(expected)) {
          ylab <- bquote(.(ylab) / .(lab.EVm))
        }
      }
    }
  }
  
  ## choose suitable ranges on the axes, unless specified by user
  if (missing(xlim)) xlim <- if (x.log) c(1, N.max) else c(0, N.max)
  if (missing(ylim)) ylim <- if (y.log) c(1, 2 * V.max) else c(0, 1.05 * V.max)

  ## set up plotting region and labels
  plot(1, 1, type="n", xlim=xlim, ylim=ylim, log=log, xaxs="i", yaxs="i",
       xlab=xlab, ylab=ylab, main=main)

  for (i in 1:n.vgc) {                # go through all specified VGCs

    curr.vgc <- VGCs[[i]]
    x.values <- N(curr.vgc)
    y.values <- if (missing(m)) V(curr.vgc) else Vm(curr.vgc, m)
    thin.lwd <- if (lwd[i] < 1) .5 else lwd[i] / 2 # for thin lines, otherwise matching style of main line

    if (x.log) x.values[x.values <= 0] <- xlim[1] / 10
    ## draw confidence intervals for expected VGC with variances 
    if (attr(curr.vgc, "hasVariances") && is.numeric(conf.level)) { 
      variance <- if (missing(m)) VV(curr.vgc) else VVm(curr.vgc, m)
      st.dev <- sqrt(variance)
      factor <- - qnorm((1 - conf.level) / 2) # approximate confidence interval
      y.minus <- y.values - factor * st.dev
      y.plus <- y.values + factor * st.dev

      if (y.log) y.minus[y.minus <= 0] <- ylim[1] / 10
      if (conf.style == "ticks") {
        .l <- length(x.values) # for "ticks" style, may need to reduce number of data points
        .idx <- 1:.l
        if (.l > 100) {
          .thin <- ceiling(.l / 100) # sample at regular intervals so that there are at most 100 data points
          .idx <- rev(seq(.l, 1, by=-.thin)) # ensure that last data point is included
        }
        segments(x.values[.idx], y.minus[.idx], x.values[.idx], y.plus[.idx], lwd=1, col=col[i])
      }
      else {
        lines(x.values, y.plus, lty=lty[i], lwd=thin.lwd, col=col[i])
        lines(x.values, y.minus, lty=lty[i], lwd=thin.lwd, col=col[i])
      }
    }

    ## draw main VGC (either V(N) or V_m(N))
      if (y.log) y.values[y.values <= 0] <- ylim[1] / 10
      lines(x.values, y.values, lty=lty[i], lwd=lwd[i], col=col[i])

    ## add VGCs for various V_m(N) (with option 'add.m')
    if (!missing(add.m)) {
      for (.m in add.m) {
        y.values <- Vm(curr.vgc, .m)
        if (y.log) y.values[y.values <= 0] <- ylim[1] / 10
        lines(x.values, y.values, lty=lty[i], lwd=thin.lwd, col=col[i])
      }
    }

  }

  if (!missing(N0)) {
    abline(v=N0, lwd=2, col="black", lty="dashed")
  }
  
  if (!missing(legend)) {             # add legend if specified by user
    legend("bottomright", inset=.02, bg="white", legend=legend, lty=lty, lwd=lwd, col=col)
  }
  
}
