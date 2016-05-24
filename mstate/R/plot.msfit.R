plot.msfit <- function(x, type=c("single", "separate"), cols,
    xlab="Time", ylab="Cumulative hazard", ylim, lwd, lty,
    legend, legend.pos, bty="n", ...)
{
  if (!inherits(x, "msfit"))
    stop("'x' must be a 'msfit' object")
  msf1 <- x$Haz
  K <- max(msf1$trans)
  msft <- unique(msf1$time) # the time points
  nt <- length(msft)
  msfp <- matrix(msf1[,2], nrow=nt) # the cumulative hazards in matrix (nt x K)
  type <- match.arg(type)
  if (missing(legend))
    legend <- to.trans2(x$trans)$transname
  if (type=="single") {
    if (missing(cols)) cols <- 1:K
    if (missing(ylim)) ylim <- c(0, max(msfp))
    if (missing(lwd)) lwd <- 1
    if (missing(lty)) lty <- rep(1, K)
    plot(msft, msfp[,1], type="s", ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1], lwd=lwd,
         lty=lty[1], ...)
    for (k in 2:K) lines(msft, msfp[,k], type="s", col=cols[k], lwd=lwd, lty=lty[k], ...)
    if (missing(legend.pos))
      legend("topleft", legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
    else
      legend(legend.pos[1], legend.pos[2], legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
  }
  else if (type=="separate") {
    if (missing(cols)) cols <- rep(1,K)
    if (missing(lwd)) lwd <- 1
    if (missing(lty)) lty <- 1
    for (k in 1:K) {
      if (missing(ylim))
        plot(msft, msfp[,k], type="s", xlab=xlab, ylab=ylab, col=cols[k], lwd=lwd, ...)
      else
        plot(msft, msfp[, k], type="s", ylim=ylim, xlab=xlab, ylab=ylab, col=cols[k],
             lwd=lwd, ...)
      title(main=legend[k])
    }
  }
  return(invisible())
}
