panel.histogram.hh <-
  if.R(r={
    function (x, breaks, equal.widths = TRUE, type = "density",
              nint = round(log2(length(x)) + 1), alpha = plot.polygon$alpha,
              col = plot.polygon$col, border = plot.polygon$border,
              lty = plot.polygon$lty, lwd = plot.polygon$lwd, ..., identifier = "histogram")
      {
        horizontal <- list(...)$horizontal
        horizontal <- if(is.null(horizontal)) FALSE else horizontal
        plot.polygon <- trellis.par.get("plot.polygon")
        xscale <- current.panel.limits()$xlim
        yscale <- current.panel.limits()$ylim  ## new
        if (!horizontal) ## default
        panel.lines(x = xscale[1] + diff(xscale) * c(0.05, 0.95),
                    y = c(0, 0), col = border, lty = lty, lwd = lwd, alpha = alpha,
                    identifier = paste(identifier, "baseline", sep = "."))
        else ## transposed
        panel.lines(y = yscale[1] + diff(yscale) * c(0.05, 0.95),
                    x = c(0, 0), col = border, lty = lty, lwd = lwd, alpha = alpha,
                    identifier = paste(identifier, "baseline", sep = "."))
        if (length(x) > 0) {
          if (is.null(breaks)) {
            breaks <- if (is.factor(x))
              `seq_len`(1 + nlevels(x)) - 0.5
            else if (equal.widths)
              do.breaks(range(x, finite = TRUE), nint)
            else quantile(x, 0:nint/nint, na.rm = TRUE)
          }
          ## h <- lattice:::hist.constructor(x, breaks = breaks, ...)
          h <- lattice.hist.constructor(x, breaks = breaks, ...)
          y <- if (type == "count")
            h$counts
          else if (type == "percent")
            100 * h$counts/length(x)
          else h$intensities
          breaks <- h$breaks
          nb <- length(breaks)
          if (length(y) != nb - 1)
            warning("problem with 'hist' computations")
          if (nb > 1) {
              if (!horizontal) ## default
                panel.rect(x = breaks[-nb], y = 0, height = y, width = diff(breaks),
                           col = col, alpha = alpha, border = border, lty = lty,
                           lwd = lwd, just = c("left", "bottom"), identifier = identifier)
              else ## transposed
                panel.rect(y = breaks[-nb], x = 0, height = diff(breaks), width = y,
                           col = col, alpha = alpha, border = border, lty = lty,
                           lwd = lwd, just = c("left", "bottom"), identifier = identifier)
          }
        }
      }
    ##<bytecode: 0x0c0f0350>
    ##<environment: namespace:lattice>
  },
       s={
         function(x, y, col = trellis.par.get("bar.fill")$col, border = 1, ...,
                  transpose=FALSE)
           {
             if (transpose) {
               z <- x
               x <- y
               y <- z
               n <- length(x)
               ymin <- rep(0, n - 1)
               ymax <- par("usr")[2]
               y <- y[-1]
               y[y > ymax] <- ymax
               polygon(c(rbind(ymin, ymin, y, y, NA)),
                       c(rbind(x[-1], x[ - n], x[ - n], x[-1], NA)),
                       col = col,
                       border = as.numeric(border), ...)
             }
             else {
               n <- length(x)
               ymin <- rep(0, n - 1)
               ymax <- par("usr")[4]
               y <- y[-1]
               y[y > ymax] <- ymax
               polygon(c(rbind(x[-1], x[ - n], x[ - n], x[-1], NA)),
                       c(rbind(ymin, ymin, y, y, NA)),
                       col = col,
                       border = as.numeric(border), ...)
             }
           }
       })
## if.R(r=assignInNamespace("panel.histogram.hh", panel.histogram.hh, "HH"),
##      s={})
## source("c:/HOME/rmh/HH-R.package/HH/R/histogram.hh.R")
## trace(panel.histogram.hh, exit=browser)


## test <- data.frame(yy=rep(1:5, 1:5))
## test.hs <- histogram(~yy, data=test, type="count")
## update(test.hs, main="vertical")

## if.R(r={
##   test.hs.t <- histogram(~yy, data=test, type="count", horizontal=TRUE, panel=panel.histogram.hh)
##   test.hs.t <- update(test.hs.t,
##                       xlim=test.hs$y.limits, ylim=test.hs$x.limits,
##                       xlab=if(test.hs$ylab==TRUE) test.hs$ylab.default else test.hs$ylab,
##                       ylab=test.hs$xlab,
##                       main="horizontal")
##   test.hs.t
## }, s={
##   t(test.hs)
## })

