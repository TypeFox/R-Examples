fillplot <- function(x,y1,y2,col,lwd)
  # y2>y1, x (ascending order), y1, y2 same length, added to existing plot, intended for type="s"
{
  nx <- length(x)
  # add mini-bit of space, this is to incorporate the possibility of a jump at the end 
  x <- c(x, x[nx]+0.1*diff(range(x))) 
  xx <- c(rep(x,c(1,rep(2,nx-1),1)),rep(rev(x),c(1,rep(2,nx-1),1)))
  yy <- c(rep(y1,rep(2,nx)),rep(rev(y2),rep(2,nx)))
  polygon(xx,yy,col=col,lwd=lwd)
}

plot.probtrans <- function(x, from=1, type=c("stacked","filled","single","separate"), ord,
                           cols, xlab="Time", ylab="Probability", xlim, ylim, lwd, lty, cex, legend, legend.pos,
                           bty="n", xaxs="i", yaxs="i", ...)
  # ord for "stacked" and "filled "only, cex for text only
{
  if (!inherits(x, "probtrans"))
    stop("'x' must be a 'probtrans' object")
  trans <- x$trans
  S <- dim(trans)[1]
  if ((from<1) | (from>S)) stop("'from' incorrect")
  pt1 <- x[[from]]
  ptt <- pt1$time # the time points
  nt <- length(ptt)
  ptp <- pt1[,2:(S+1)] # those are the actual transition probabilities
  if (missing(ord)) ord <- 1:S
  if (missing(lwd)) lwd <- 1
  lwd <- rep(lwd, S)
  lwd <- lwd[1:S]
  if (missing(lty)) lty <- 1
  lty <- rep(lty, S)
  lty <- lty[1:S]
  if (missing(cex)) cex <- 1
  cex <- rep(cex, S)
  cex <- cex[1:S]
  type <- match.arg(type)
  if (missing(cols)) {
    if (type=="filled" | type=="single") cols <- rev(RColorBrewer::brewer.pal(S, "RdYlGn"))
    else cols <- 1
  }
  cols <- rep(cols, S)
  cols <- cols[1:S]
  if (missing(legend)) {
    legend <- dimnames(trans)[[2]]
    if (is.null(legend)) legend <- as.character(1:S)
  }
  else if (length(legend) != S) stop("legend has incorrect length")
  if (type=="single") {
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,max(ptp))
    plot(ptt, ptp[,1], type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], lty=lty[1], ...)
    for (s in 2:S) lines(ptt, ptp[,s], type="s", col=cols[s], lwd=lwd[s], lty=lty[s], ...)
    if (missing(legend.pos))
      legend("topright", legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
    else
      legend(legend.pos[1], legend.pos[2], legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
  }
  else if (type=="stacked") {
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,1)
    # finally change order according to ord
    lwd <- lwd[ord]; lty <- lty[ord]; cex <- cex[ord]; cols <- cols[ord]
    y0 <- 0
    ptpsum <- ptp[,ord[1]]
    eps <- (xlim[2]-xlim[1])/50
    nt <- sum(ptt <= xlim[2]-eps)
    dy <- ptp[nt,ord[1]]
    y <- y0 + dy/2
    y1 <- y0 + dy
    plot(ptt, ptpsum, type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], ...)
    text(xlim[2]-eps, y, legend[ord[1]], adj=1, cex=cex)
    for (s in 2:S) {
      ptpsum <- ptpsum + ptp[,ord[s]]
      lines(ptt, ptpsum, type="s", col=cols[s], lwd=lwd[s], ...)
      y0 <- y1
      dy <- ptp[nt,ord[s]]
      y <- y0 + dy/2
      y1 <- y0 + dy
      text(xlim[2]-eps, y, legend[ord[s]], adj=1, cex=cex)
    }
  }
  else if (type=="filled") {
    par(xaxs=xaxs, yaxs=yaxs)
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,1)
    # finally change order according to ord
    lwd <- lwd[ord]; lty <- lty[ord]; cex <- cex[ord]; cols <- cols[ord]
    y0 <- 0
    eps <- (xlim[2]-xlim[1])/50
    nt <- sum(ptt <= xlim[2])
    ptt <- ptt[1:nt]
    ptplow <- rep(0,nt)
    ptpup <- ptp[1:nt, ord[1]]
    dy <- ptp[nt, ord[1]]
    y <- y0 + 0.5*dy
    y1 <- y0 + dy
    plot(ptt, ptpup, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], ...)
    fillplot(ptt, ptplow, ptpup, col=cols[1],lwd=lwd[1])
    text(xlim[2]-eps, y, legend[ord[1]], adj=1, cex=cex[1])
    for (s in 2:S) {
      ptplow <- ptpup
      ptpup <- ptpup + ptp[1:nt,ord[s]]
      fillplot(ptt, ptplow, ptpup, col=cols[s],lwd=lwd[s])
      y0 <- y1
      dy <- ptp[nt, ord[s]]
      y <- y0 + 0.5*dy
      y1 <- y0 + dy
      text(xlim[2]-eps, y, legend[ord[s]], adj=1, cex=cex[s])
    }
    box()
  }
  else if (type=="separate") {
    if (missing(xlim)) xlim <- range(ptt)
    for (s in 1:S) {
      if (missing(ylim)) # each on its own y-scale
        plot(ptt, ptp[,s], type="s", xlim=xlim,  xlab=xlab, ylab=ylab, col=cols[s], lwd=lwd[s])
      else
        plot(ptt, ptp[,s], type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[s],
             lwd=lwd[s], ...)
      title(main=legend[s])
    }
  }
  return(invisible())
}
