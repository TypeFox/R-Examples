plot.tfl <- function (x, y, ...,
                      min.rank=1, max.rank=NA, log=c("","x","y","xy"), 
                      type=c("p","l","b","c","o","h","s"),
                      bw=zipfR.par("bw"), cex=1,
                      xlim=NULL, ylim=NULL, 
                      xlab="rank", ylab="frequency", legend=NULL,
                      main="Type-Frequency List (Zipf ranking)",
                      pch=NULL, lty=NULL, lwd=NULL, col=NULL)
{
  ## collect all specified TFL in single list
  TFLs <- list(x)   # this is a bit complicated because of the plot() prototype
  if (! missing(y)) {
    TFLs <- c(TFLs, list(y), list(...))
  }
  n.tfl <- length(TFLs)
  
  ## check other arguments
  if (length(log) > 1 || log != "") log <- match.arg(log)
  type <- match.arg(type)
  if (!missing(legend) && length(legend) != n.tfl) stop("'legend' argument must be character or expression vector of same length as number of VGCs")
  if (any(sapply(TFLs, function (.TFL) attr(.TFL, "incomplete")))) stop("plotting of incomplete type-frequency lists is not supported")
  x.log <- log == "x" || log == "xy"
  y.log <- log == "y" || log == "xy"

  ## determine range of ranks available & check requested min/max
  max.rank.available <- max(sapply(TFLs, function (.TFL) length(.TFL$f)))
  if (missing(max.rank)) {
    max.rank <- max.rank.available
  } else {
    if (min.rank >= max.rank.available) stop(sprintf("no data available for ranks > %d (min.rank)", min.rank))
  }
  if (min.rank < 1) stop(sprintf("min.rank=%d invalid, must be integer >= 1", min.rank))
  if (min.rank >= max.rank) stop(sprintf("min.rank=%d must be smaller than max.rank=%d", min.rank, max.rank))
  ranks <- min.rank:max.rank

  ## collect list of Zipf-ranked type frequency vectors, then determine default range for y-axis
  freqs <- lapply(TFLs, function (.TFL) sort(.TFL$f, decreasing=TRUE))
  f.max <- max(unlist(lapply(freqs, function (.F) na.omit(.F[ranks]))))
  
  ## get default styles unless manually overridden
  if (missing(pch)) pch <- zipfR.par("pch", bw.mode=bw)
  if (missing(lty)) lty <- zipfR.par("lty", bw.mode=bw)
  if (missing(lwd)) lwd <- zipfR.par("lwd", bw.mode=bw)
  if (missing(col)) col <- zipfR.par("col", bw.mode=bw)
  
  ## choose suitable ranges on the axes, unless specified by user
  if (missing(xlim)) xlim <- c(min.rank, max.rank)
  if (missing(ylim)) ylim <- if (y.log) c(2/3, 1.5*f.max) else c(0, 1.05 * f.max)

  ## set up plotting region and labels
  plot(1, 1, type="n", xlim=xlim, ylim=ylim, log=log, xaxs="r", yaxs="i",
       xlab=xlab, ylab=ylab, main=main)

  for (i in 1:n.tfl) {                # go through all specified TFLs
    curr.freqs <- freqs[[i]]
    curr.ranks <- ranks[!is.na(curr.freqs[ranks])] # reduce ranks to valid set
    if (length(curr.ranks) >= 1) {
      if (type == "s") {
        idx <- c(TRUE, diff(curr.freqs[curr.ranks]) != 0)
        idx[length(idx)] <- TRUE # extend last horizontal line to end of data
        points(curr.ranks[idx], curr.freqs[curr.ranks[idx]], type="s", cex=cex, pch=pch[i], lty=lty[i], lwd=lwd[i], col=col[i])
      } else {
        points(curr.ranks, curr.freqs[curr.ranks], type=type, cex=cex, pch=pch[i], lty=lty[i], lwd=lwd[i], col=col[i])
      }
    }
  }

  if (!missing(legend)) {             # add legend if specified by user
    legend.args <- list("topright", inset=.02, bg="white", legend=legend, col=col)
    if (type %in% c("p","b","o")) legend.args <- append(legend.args, list(pch=pch, pt.cex=1.4*cex, pt.lwd=lwd))
    if (type != "p") legend.args <- append(legend.args, list(lwd=lwd))
    do.call("legend", legend.args)
  }  
}
