

# x: results from aggregateColor()
aggregateColorPlot <- function(x, print.label=TRUE, label.font=1, label.cex=0.65, buffer.pct=0.02, print.n.hz=FALSE, rect.border='black', horizontal.borders=FALSE, ...) {
 
  # extract just the scaled data from the results of aggregateColor()
  s.scaled <- x$scaled.data
  
  # get max re-scaled summation for xlim
  max.plot <- max(sapply(s.scaled, function(i) sum(i$weight)))
  
  # setup plot
  plot(1,1, type='n', xlim=c(0, max.plot), ylim=c(length(names(s.scaled))+0.5, 0.5), axes=FALSE, ylab='', xlab='Cumulative Proportion', ...)
  # iterate over horizons
  for(i in seq_along(names(s.scaled))) {
    s.i <- s.scaled[[i]]
    n.colors <- nrow(s.i)
    if(n.colors > 0) {
      # get an index to the last weight
      last.weight <- length(s.i$weight)
      # compute cumulative left / right rectangle boundaries
      x.left <- cumsum(c(0, s.i$weight[-last.weight]))
      x.right <- c(x.left[-1], x.left[last.weight] + s.i$weight[last.weight])
      
      # plot rectanges from vectorized coordinates / colors
      rect(xleft=x.left, ybottom=i-0.5, xright=x.right, ytop=i+0.5, col=s.i$soil_color, border=rect.border)
      
      # compute center point for color labels
      centers <- (x.right + x.left) / 2
      
      # create label
      if(print.n.hz)
        color.labels <- paste0(s.i$munsell, '\n', '(', s.i$n.hz, ')')
      else
        color.labels <- s.i$munsell
      
      # determine if there is enough room to label colors: some % of visible space on plot
      # first, get the plot aspect ratio
      plot.w <- par("pin")[1]/diff(par("usr")[1:2])
      plot.h <- par("pin")[2]/diff(par("usr")[3:4])
      plot.asp <- abs(plot.w / plot.h)
      # get text heigths, as we will be printing labels at 90 deg
      text.heights <- abs(strheight(color.labels, cex=label.cex, font=label.font))
      # convert text heights into equivelent widths
      text.heights <- text.heights / plot.asp
      # compare re-scaled text heights with rectangle widths (weights) + some buffer
      label.fits <- which(text.heights < (s.i$weight - buffer.pct) )
      
      # print labels
      if(print.label & (length(label.fits) > 0))
        text(x=centers[label.fits], y=i, labels=color.labels[label.fits], col='white', font=label.font, cex=label.cex, srt=90)
    }
  }
  
  # add horizontal separator lines, typically used when rectange borders are not drawn
  if(horizontal.borders){
    hz.line.y <- 1:(length(names(s.scaled))-1) + 0.5
    segments(x0 = 0, y0 = hz.line.y, x1 = 1, y1 = hz.line.y, lwd=2)
  }
  
  ## TODO: make these adjustable
  # add axes
  axis(1, at=round(seq(0, 1, length.out = 11), 2))
  axis(2, at = seq_along(names(s.scaled)), labels = names(s.scaled), las=2, tick=FALSE, font=2, hadj=1, line=-2.125, cex.axis=1)
  
}
