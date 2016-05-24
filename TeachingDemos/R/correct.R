cor.rect.plot <- function(x,y, corr=TRUE,
                    xlab=deparse(substitute(x)),
                    ylab=deparse(substitute(y)),
                    col=c('#ff000055','#0000ff55'),
                    ... ) {
  xy <- xy.coords(x,y,xlab=xlab, ylab=ylab)
	xm <- mean(xy$x)
	ym <- mean(xy$y)

  op <- par(mar=c(5,4,4,4))
  on.exit(par(op))

  plot(xy$x, xy$y, xlab=xlab, ylab=ylab, pty='s')

  xt <- scale(xy$x, scale=corr)
  yt <- scale(xy$y, scale=corr)

  xtt <- pretty(xt)
  ytt <- pretty(yt)

  xut <- if(corr) {
    xtt * attr(xt, 'scaled:scale') + attr(xt, 'scaled:center')
  } else {
    xtt + attr(xt, 'scaled:center')
  }
  yut <- if(corr) {
    ytt * attr(yt, 'scaled:scale') + attr(yt, 'scaled:center')
  } else {
    ytt + attr(yt, 'scaled:center')
  }

  axis(3,at=xut, labels=xtt)
  axis(4,at=yut, labels=ytt)

  abline(h=ym)
  abline(v=xm)


  ord <- order( xt^2+yt^2, decreasing=TRUE )
  w <- xt[ord]*yt[ord] > 0
  rect(xm, ym, xy$x[ord][w], xy$y[ord][w], col= col[1] )
  rect(xm, ym, xy$x[ord][!w], xy$y[ord][!w], col= col[2] )
  points(xy$x,xy$y)
}

