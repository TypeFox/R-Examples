mdcplot <- function(x, rescale = TRUE, ...) {

  ## check if x is a single data.frame as returned by the old version
  ## of mdcc, or a list of data.frames as returned by the newer
  ## version (> 1.2.2); if the first is true, only take the coefs
  if (!is.data.frame(x)) {
    x <- x$coef
  }

  blues <- colorRamp(c("#FFFFFF", "#395cd4"))
  reds <- colorRamp(c("#FFFFFF", "#dd291c"))

  m <- dim(x)[1]                        # parameters
  n <- dim(x)[2]                        # time steps

  xs <- x[m:1,]                         # invert rows for plotting "bottom-up"
    
  if (rescale) {
    pos.max <- max(x)
    neg.max <- abs(min(x))
  }
  
  ## plot frame

  op <- par(no.readonly = TRUE)
  par(oma = c(2,3,0.1,5), mai = c(1, 0.8, 0.2, 1))
  plot(c(0.4, n + 0.4), c(0.4, m + 0.4), type = "n", xaxt = "n", yaxt
       = "n", ylab = "", xlab = "", ...)

  axis(side = 1, at = 1:n, labels = colnames(x), las = 2)
  axis(side = 2, at = 1:m, labels = rev(rownames(x)), las = 1)

  for (i in 1:m) {
    for (j in 1:n) {
      if (xs[i,j] <= 0) {
        
        polygon(c(j-0.4, j+0.4, j+0.4, j-0.4), c(i-0.4, i-0.4, i+0.4,
                                                 i+0.4), lty = 0, col =
                ifelse(rescale,
                       rgb(reds(abs(xs[i,j])/neg.max),
                           maxColorValue = 255),
                       rgb(reds(abs(xs[i,j])),
                           maxColorValue = 255)))

      } else {

        polygon(c(j-0.4, j+0.4, j+0.4, j-0.4), c(i-0.4, i-0.4, i+0.4,
                                                 i+0.4), lty = 0, col =
                ifelse(rescale,
                       rgb(blues(xs[i,j]/pos.max),
                           maxColorValue = 255),
                       rgb(blues(xs[i,j]),
                           maxColorValue = 255)))

      }
    }
  }

  par(xpd = NA)

  ## draw color legend

  ## divide total span of parameter instances into 11 parts (5
  ## negative, 5 positive, 1 zero)

  leg.unit <- (m/11)
  start.unit <- leg.unit/4 - leg.unit
  right.pos <- n + n/10
  leg.width <- n/15
  values <- seq(-1, 1, length = 11)
  neg.rescaled.values <- round(seq(min(x), 0, length = 6), 2)
  pos.rescaled.values <- rev(round(seq(max(x), 0, length = 6), 2)[-6])
  rescaled.values <- c(neg.rescaled.values, pos.rescaled.values)
  
  for (i in 1:11) {
    if (values[i] <= 0) {
      polygon(c(right.pos, right.pos+leg.width, right.pos+leg.width, right.pos),
              c(start.unit + ((i-1)*leg.unit), start.unit + ((i-1)*leg.unit),
                start.unit + (i*leg.unit), start.unit +
                (i*leg.unit)), col = rgb(reds(abs(values[i])),
                                 maxColorValue = 255), lty = 0)
      text(right.pos + leg.width + 1, start.unit + (i*leg.unit) -
           leg.unit/2, ifelse(rescale, rescaled.values[i], values[i]),
           pos = 4)
    } else {
      polygon(c(right.pos, right.pos+leg.width, right.pos+leg.width, right.pos),
              c(start.unit + ((i-1)*leg.unit), start.unit + ((i-1)*leg.unit),
                start.unit + (i*leg.unit), start.unit +
                (i*leg.unit)), col = rgb(blues(values[i]),
                                 maxColorValue = 255), lty = 0)
      text(right.pos + leg.width + 1, start.unit + (i*leg.unit) -
           leg.unit/2, ifelse(rescale, rescaled.values[i], values[i]),
           pos = 4)
    }
  }
  par(op)
}
