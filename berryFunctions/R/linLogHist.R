#' lin-log transition histogram
#' 
#' Draw histograms that gradually transform from a linear to a logarithmic axis (animation)
#' 
#' @return Returned invisibly: transformation values used. Plotted: \code{steps} number of images.
#' @note It's best to save the plots into a pdf or wrap it within\cr
#'       \code{png("Transition\%03d"); linLogHist(x); dev.off()}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, April 2015
#' @seealso \code{\link{linLogTrans}}
#' @keywords dplot hplot dynamic
#' @export
#' @examples
#' 
#' x <- rlnorm(700, m=3)
#' hist(x, col=4)
#' hist(log10(x), xaxt="n"); logAxis(1); hist(log10(x), col=4, add=TRUE)
#' 
#' op <- par()
#' linLogHist(x, steps=8, sleep=0.01) # 0.05 might be smoother
#' 
#' linLogHist(x, xlab="ddd", breaks=30, steps=3, write_t=FALSE, yaxt="n", freq=FALSE,
#'    main="", parexpr='par(mar=c(2,0.5,1.5,0.5), mgp=c(1.8,1,0))',
#'    endexpr='mtext("Probability Density", line=-1.2, adj=0.03, outer=T)')
#' par(op)
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices such as pdf,
#' ## so this example is excluded from running in the checks.
#' pdf("LinLogTransitionAnimation.pdf")
#' linLogHist(x, main="Example Transition", steps=20, freq=FALSE)
#' dev.off()
#' 
#' # if you have FFmpeg installed, you can use the animation package like this:
#' library2(animation)
#' saveVideo(linLogHist(x, steps=50), video.name="linlog_anim.mp4", interval=0.08,
#' ffmpeg="C:/ffmpeg-20150424-git-cd69c0e-win64-static/bin/ffmpeg.exe")
#' }
#' 
#' @param x x values to be plotted in animation
#' @param steps Number of steps in transition. DEFAULT: 100
#' @param breaks \code{\link{hist}} breaks. DEFAULT: 20
#' @param col \code{\link{hist}} color. DEFAULT: "blue"
#' @param las \code{\link{par}} LabelAxisStyle (numbers upright). DEFAULT: 1
#' @param xlab Label for the x axis. DEFAULT: deparse(substitute(x))
#' @param xlim xlim range in non-log units. DEFAULT: range(x, finite=TRUE)
#' @param box Draw box at the end to overplot \code{\link{abline}s} crossing the box? DEFAULT: TRUE
#' @param parexpr Characterized Expression to set \code{\link{par}}, eg. \code{parexpr='par(mar=c(2,0.5,1.5,0.5), mpg=c(1.8,1,0))'}
#' @param endexpr Characterized Expression executed at the end of the plot, eg. \code{endexpr='mtext("Probability Density", line=-1, adj=0.03, outer=T)'}
#' @param axisargs List of arguments passed to \code{\link{logVals}}, like base. DEFAULT: NULL
#' @param sleep Pause time between frames, in seconds, passed to \code{\link{Sys.sleep}}. DEFAULT: 0
#' @param axisargs2 List of arguments passed to \code{\link{logAxis}} in the final plot. DEFAULT: NULL
#' @param firstplot plot on linear scale first? DEFAULT: TRUE
#' @param lastplot plot on logarithmic scale at the end? DEFAULT: TRUE
#' @param write_t write transformation value in lower right corner? DEFAULT: TRUE
#' @param values_t Supply vector with values for transformation (1/t). Overides steps. 
#'        If you have a better algorithm than I do, please let me know! DEFAULT: NULL
#' @param \dots further arguments passed to \code{\link{hist}}, like freq, main, xlim, ylab. Excluded: x, xaxt, possibly add
#' 
linLogHist <- function(
x,
steps=100,
breaks=20,
col="blue",
las=1,
xlab=deparse(substitute(x)),
xlim=range(x, finite=TRUE),
box=TRUE,
parexpr,
endexpr,
sleep=0,
axisargs=NULL,
axisargs2=NULL,
firstplot=TRUE,
lastplot=TRUE,
write_t=TRUE,
values_t=NULL,
...)
{
# x must be deparsed before it's evaluated (or something like that)
xlab <- xlab
# Tansformation values ---------------------------------------------------------
allt <- if(is.null(values_t))  linLogTrans(x, steps=steps, plot=FALSE)  else  values_t
# Plot the histograms ----------------------------------------------------------
# Plot on linear scale first:
if(firstplot)
  {
  if(!missing(parexpr)) eval(parse(text=parexpr))
  hist(x, breaks=breaks, col=col, las=las, xlab=xlab, xlim=xlim, ...)
  if(box) graphics::box("plot")
  if(!missing(endexpr)) eval(parse(text=endexpr))
  }
#
# Log labels and lines:
lv <- do.call(logVals, owa(list(from=x), axisargs))
# Images:
for(t in allt)
  {
  # Plot single frame:
  if(!missing(parexpr)) eval(parse(text=parexpr))
  hist(x^(1/t), breaks=breaks, col=col, las=las, xlab=xlab, xaxt="n", xlim=xlim^(1/t), ...)
  # draw grey lines at 10^n values and label appropriate ones:
  abline(v=(lv$all)^(1/t), col=8)
  axis(1, (lv$vals)^(1/t), lv$labs, las=las)
  hist(x^(1/t), breaks=breaks, col=col, add=TRUE, ...)
  if(box) graphics::box("plot")
  # write transformation value:
  if(write_t) title(sub=paste("t =", sprintf("%6.2f", t)), adj=1)
  if(!missing(endexpr)) eval(parse(text=endexpr))
  # slow down frame passing:
  if(sleep!=0) Sys.sleep(sleep)
  } # End for loop
# Final image
if(lastplot)
  {
  if(!missing(parexpr)) eval(parse(text=parexpr))
  hist(log10(x), breaks=breaks, col=col, las=las, xlab=xlab, xaxt="n", xlim=log10(xlim), ...)
  do.call(logAxis, args=owa(c(side=1, box=box), axisargs2))
  hist(log10(x), breaks=breaks, col=col, add=TRUE, ...)
  if(!missing(endexpr)) eval(parse(text=endexpr))
  }
return(invisible(allt))
} # end of function ------------------------------------------------------------
