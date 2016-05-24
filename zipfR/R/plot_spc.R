plot.spc <- function(x, y, ...,                                           
                     m.max=if (log=="") 15 else 50,
                     log="", conf.level=.95,
                     bw=zipfR.par("bw"), points=TRUE,
                     xlim=NULL, ylim=NULL,
                     xlab="m", ylab="V_m", legend=NULL,
                     main="Frequency Spectrum",
                     barcol=NULL, pch=NULL, lty=NULL, lwd=NULL, col=NULL)
{
  ## collect all specified frequency spectra in single list
  spectra <- list(x)   # this is a bit complicated because of the plot() prototype
  if (! missing(y)) {
    spectra <- c(spectra, list(y), list(...))
  }
  n.spc <- length(spectra)
  
  ## check other arguments & collect some statistics 
  if (! (log %in% c("", "x", "y", "xy"))) stop("allowed values for 'log' argument are '', 'x', 'y' and 'xy'")
  V.max <- max(sapply(spectra, function (.S) max(Vm(.S, 1:m.max))))
  if (!missing(legend) && length(legend) != n.spc)
    stop("'legend' argument must be character or expression vector of same length as number of spectra")
  x.log <- log == "x" || log == "xy"
  y.log <- log == "y" || log == "xy"
  
  ## get default styles unless manually overridden
  if (missing(barcol)) barcol <- zipfR.par("barcol", bw.mode=bw)
  if (missing(pch)) pch <- zipfR.par("pch", bw.mode=bw)
  if (missing(lty)) lty <- zipfR.par("lty", bw.mode=bw)
  if (missing(lwd)) lwd <- zipfR.par("lwd", bw.mode=bw)
  if (missing(col)) col <- zipfR.par("col", bw.mode=bw)

  ## typeset default label on y-axis, depending on whether spectra are observed or expected
  expected <- sapply(spectra, function (.S) attr(.S, "expected"))
  if (missing(ylab)) {
    if (all(expected)) {
      ylab <- expression(E * group("[", V[m], "]"))
    }
    else if (all(!expected)) {
      ylab <- expression(V[m])
    }
    else {
      ylab <- expression(V[m] / E * group("[", V[m], "]"))
    }
  }

  ## choose suitable ranges on the axes, unless specified by user
  if (missing(xlim)) xlim <- c(1, m.max)
  if (missing(ylim)) ylim <- if (y.log) c(0.1, 2 * V.max) else c(0, 1.05 * V.max)

  ## default: non-logarithmic barplot (log parameter not specified)
  if (missing(log)) {
    my.data <- sapply(spectra, function (.S) Vm(.S, 1:m.max))
    barplot(t(my.data), beside=TRUE, ylim=ylim,
            col=barcol[1:n.spc], names.arg = 1:m.max,
            xlab=xlab, ylab=ylab, main=main, legend=legend)
  }
  ## lines or points+lines for any kind of logarithmic plot (log="" for non-logarithmic scale)
  else {
    plot(1, 1, type="n", xlim=xlim, ylim=ylim, log=log, yaxs="i",
         xlab=xlab, ylab=ylab, main=main)

    for (i in 1:n.spc) {                # go through all specified frequency spectra
      curr.spc <- spectra[[i]]
      x.values <- 1:m.max
      y.values <- Vm(curr.spc, 1:m.max)

      ## draw confidence intervals for expected spectrum with variances 
      if (attr(curr.spc, "hasVariances") && is.numeric(conf.level)) { 
        st.dev <- sqrt(VVm(curr.spc, 1:m.max))
        factor <- - qnorm((1 - conf.level) / 2) # approximate confidence interval
        y.minus <- y.values - factor * st.dev
        y.plus <- y.values + factor * st.dev

        if (y.log) y.minus[y.minus <= 0] <- ylim[1] / 10
        segments(x.values, y.minus, x.values, y.plus, lwd=1, col=col[i])
        segments(x.values-.4, y.plus, x.values+.4, y.plus, lwd=1, col[i])
        segments(x.values-.4, y.minus, x.values+.4, y.minus, lwd=1, col[i])
      }

      if (y.log) y.values[y.values <= 0] <- ylim[1] / 10
      if (points) {                     # overplot points and thin lines
        points(x.values, y.values, type="o",
               pch=pch[i], lty="solid", col=col[i])
      }
      else {                            # lines only, with different styles
        lines(x.values, y.values, lty=lty[i], lwd=lwd[i], col=col[i])
      }
    }

    if (!missing(legend)) {             # add legend if specified by user
      if (points) {
        legend("topright", inset=.02, bg="white", legend=legend, pch=pch, col=col)
      }
      else {
        legend("topright", inset=.02, bg="white", legend=legend, lty=lty, lwd=lwd, col=col)
      }
    }
  }
  
}
