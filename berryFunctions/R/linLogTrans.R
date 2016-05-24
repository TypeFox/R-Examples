#' Animation for transition from linear to logarithmic axis
#' 
#' draw images that gradually transform from a linear to a logarithmic axis
#' 
#' @return Returned invisibly: transformation values used. Plotted: \code{steps} number of images.
#' @note if(steps>1000) steps <- 1000. In the unlikely case you need more steps, please let me know and I'll change the code.\cr 
#'       It's best to save the plots into a pdf (see the example) or wrap it within\cr
#'       \code{png("Transition\%03d"); linLogTrans(x,y); dev.off()}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, June 2014
#' @seealso \code{\link{logVals}}
#' @references x^(1/t) is based on the first comment on \url{http://stackoverflow.com/questions/15994442/}\cr 
#'   besides the nice graphic properties of logtransformations, check this page for the implications on rates of change: \cr
#'   \url{http://sfew.websitetoolbox.com/post/show_single_post?pid=1282690259&postcount=4}\cr
#'   \url{http://sfew.websitetoolbox.com/post/show_single_post?pid=1282691799&postcount=5}
#' @keywords dplot hplot dynamic
#' @export
#' @examples
#'  
#' set.seed(42);  x <- 10^rnorm(100, 3);  y <- runif(100)
#' linLogTrans(x,y, steps=15, sleep=0.01) # 0.05 might be smoother...
#' linLogTrans(x,y, steps=15, log="y", ylim=c(0.1, 0.8), base=c(1,2,5))
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices such as pdf,
#' ## so this example is excluded from running in the checks.
#' pdf("LinLogTransitionAnimation.pdf")
#' linLogTrans(x,y, main="Example Transition")
#' dev.off()
#' 
#' # if you have FFmpeg installed, you can use the animation package like this:
#' library2(animation)
#' saveVideo(linLogTrans(x,y, steps=300), video.name="linlog_anim.mp4", interval=0.01,
#'     ffmpeg="C:/ffmpeg-20150424-git-cd69c0e-win64-static/bin/ffmpeg.exe")
#' 
#' 
#' # old t values were dependent on the value of steps
#' findt <- function(steps) {
#'   # t-values for x^(1/t):
#'   allt <- 10^(seq(0,2.5,len=1e4) )
#'   # selection at upper half of these values;
#'   # Otherwise, the animation slows down too much at the end
#'   f <- 1.4 # multiplication factor due to length loss by unique
#'   sel <- round(seq(1, 10, len=f*steps)^4)   #0.5*seq(1, 100, len=1.3*steps)^2 + 0.5*
#'   sel2 <- unique(round(log10(seq(1, 10, len=f*steps))*f*steps))
#'   sel2[1] <- 1
#'   sel <- sel[sel2]
#'   # final t-values for transition:
#'   allt <- unique(round(allt[sel], 2))
#'   data.frame(x=seq(1,1000,len=length(allt)), t=allt)
#'   }
#' 
#' plot(findt(1000), type="l", log="y", las=1)
#' for(i in 5:999) lines(findt(i), col=rainbow2(1000)[i])
#' d <- findt(300)
#' lines(d) # good average
#' 
#' plot(d$x[-1], diff(d$t), type="l", ylim=c(3e-3,10), yaxt="n", log="y", main="t value growth rate")
#' logAxis(2) ; lines(d$x[-1], diff(d$t))
#' d2 <- findt(1000)
#' lines(d2$x[-1], diff(d2$t), col=2)
#' lines(2:1000, diff(linLogTrans(1,1, steps=1000, plot=F)), col=4)
#' 
#' 
#' d <- findt(300)
#' pdf("degreepoly.pdf")
#' for(i in 5:30)
#'   {
#'   plot(d, log="y", type="l", lwd=3, main=i, xlim=c(0,300), ylim=c(1,2))
#'   modell <- lm(t ~  poly(x,i, raw=T), data=d)
#'   lines(x2, predict(modell, data.frame(x=1:1300)), col=2)
#'   }
#' dev.off()   # 17 is good
#' 
#' cf <- coef(lm(t ~  poly(x,17, raw=T), data=d)) # these are currently used in the function
#' x <- 1:1000
#' y <- rowSums(sapply(1:18, function(i) cf[i]*x^(i-1)), na.rm=TRUE)
#' lines(x, y, lwd=3)
#' y[1] <- 1
#' plot(x, round(y, 3), ylim=c(1,3), xlim=c(0,500), type="l", log="")
#' dput(round(y, 3))
#' 
#' findn <- function(steps) nrow(findt(steps))
#' plot(1:1000, sapply(1:1000, findn), type="l")
#' abline(b=1, a=0)
#' 
#' }
#' 
#' @param x x values to be plotted in animation
#' @param y Vector with corresponding y values
#' @param log Which axis is logarithmic, "x" or "y". DEFAULT: "x"
#' @param steps Number of steps (images) in transition (About 30\% are taken out). DEFAULT: 100
#' @param base Base passed to \code{\link{logVals}}. DEFAULT: 1
#' @param las \code{\link{par}} LabelAxisStyle (numbers upright). DEFAULT: 1
#' @param plot Plot animations at all? False to just get the t-vector (used in \code{\link{linLogHist}}). DEFAULT: TRUE
#' @param xlim xlim range in non-log units. DEFAULT: range(x, finite=TRUE)
#' @param ylim ylim range in non-log units. DEFAULT: range(y, finite=TRUE)
#' @param box Draw box at the end to overplot \code{\link{abline}s} crossing the box? DEFAULT: TRUE
#' @param parexpr Characterized Expression to set \code{\link{par}}, eg. \code{parexpr='par(mar=c(2,0.5,1.5,0.5), mpg=c(1.8,1,0))'}
#' @param endexpr Characterized Expression executed at the end of the plot, eg.
#'        \code{endexpr='mtext("Probability density", line=-1, adj=0.03, outer=T)'}
#' @param sleep Pause time between frames, in seconds, passed to \code{\link{Sys.sleep}}. DEFAULT: 0
#' @param firstplot Plot data on linear axis as additional first image? DEFAULT: TRUE
#' @param lastplot Plot data on logarithmic axis as additional last image? DEFAULT: TRUE
#' @param write_t Write transformation value in lower right corner? DEFAULT: TRUE
#' @param values_t Supply vector with values for transformation (1/t). Overides steps. If you have a better algorithm than I do, please let me know! DEFAULT: NULL for internal calculation based on size of steps.
#' @param pointsarg List of further arguments passed to points, like pch, cex, col. DEFAULT: NULL
#' @param \dots Further arguments passed only to plot, like main, xlim, ylab. Excluded: x, y, las, xaxt, type
#' 
linLogTrans <- function(
x,
y,
log="x",
steps=100,
base=1,
las=1,
plot=TRUE,
xlim=range(x, finite=TRUE),
ylim=range(y, finite=TRUE),
box=TRUE,
parexpr,
endexpr,
sleep=0,
firstplot=TRUE,
lastplot=TRUE,
write_t=TRUE,
values_t=NULL,
pointsarg=NULL,
...)
{
# Tansformation values ---------------------------------------------------------
if(is.null(values_t)) # if it's not given by user, use internal calculation:
  {
  if(length(steps)>1) steps <- max(steps, na.rm=TRUE)
  if(steps>1000) steps <- 1000
  # t-values for x^(1/t):
  # coefficients of degree 17 polynomial, see examples
  tvc <- c(0.995726006684206, 0.00310820777426466, -2.7073341927363e-05,
  6.11849831959088e-07, -7.05912829318337e-09, 4.82269115381641e-11,
  -2.02029859969675e-13, 5.30790764027315e-16, -8.53304767303152e-19,
  7.29652239791065e-22, -8.04764444262045e-26, -4.35519950517021e-28,
  3.26048883565918e-31, NA, -6.70898382748396e-38, NA, 1.04885136585542e-44,NA)
  tvx <- 1:1000 # t x values
  tvy <- rowSums(sapply(1:18, function(i) tvc[i]*tvx^(i-1)), na.rm=TRUE)
  tvy[1] <- 1
  # final t-values for transition:
  allt <- tvy[seq(1,1000, length=steps)]
  }
  else allt <- values_t
# return t values only:
if(!plot) return(allt)
# Plot the images --------------------------------------------------------------
# Plot on linear scale first:
if(firstplot)
  {
  if(!missing(parexpr)) eval(parse(text=parexpr))
  plot(x, y, las=las, type="n", xlim=xlim, ylim=ylim, ...)
  do.call(points, args=owa(d=list(x=x, y=y), a=pointsarg, "x", "y"))
  if(box) graphics::box("plot")
  if(!missing(endexpr)) eval(parse(text=endexpr))
  }
# in case people capitalize log:
log <- tolower(log)
if(log=="x") # -----------------------------------------------------------------
  {
  # Log labels and lines:
  lv <- logVals(x, base=base)
  # Images:
  for(t in allt)
     {
     if(!missing(parexpr)) eval(parse(text=parexpr))
     # Plot single frame:
     plot(x^(1/t), y, las=las, xaxt="n", type="n", xlim=xlim^(1/t), ylim=ylim, ...)
     # draw grey lines at 10^n values and label appropriate ones:
     abline(v=(lv$all)^(1/t), col=8)
     axis(1, (lv$vals)^(1/t), lv$labs, las=las)
     # user-specified arguments for points:
     pargs <- owa(d=list(x=x^(1/t), y=y), a=pointsarg, "x", "y")
     # draw original points:
     do.call(points, args=pargs)
     if(box) graphics::box("plot")
     if(!missing(endexpr)) eval(parse(text=endexpr))
     # write transformation value:
     if(write_t) title(sub=paste("t =", sprintf("%6.2f", t)), adj=1)
     # slow down frame passing:
     if(sleep!=0) Sys.sleep(sleep)
     } # End for loop
  # Final image
  if(lastplot)
    {
    if(!missing(parexpr)) eval(parse(text=parexpr))
    plot(x, y, las=las, xaxt="n", type="n", log="x", xlim=xlim, ylim=ylim, ...)
    abline(v=lv$all, col=8)
    axis(1, lv$vals, lv$labs, las=las)
    do.call(points, args=owa(d=list(x=x, y=y), a=pointsarg, "x", "y"))
    if(box) graphics::box("plot")
    if(!missing(endexpr)) eval(parse(text=endexpr))
    }
  }
else if(log=="y") # ------------------------------------------------------------
  {
  lv <- logVals(y, base=base)
  for(t in allt)
     {
     if(!missing(parexpr)) eval(parse(text=parexpr))
     plot(x, y^(1/t), las=las, yaxt="n", type="n", xlim=xlim, ylim=ylim^(1/t), ...)
     abline(h=(lv$all)^(1/t), col=8)
     axis(2, (lv$vals)^(1/t), lv$labs, las=las)
     pargs <- owa(d=list(x=x, y=y^(1/t)), a=pointsarg, "x", "y")
     do.call(points, args=pargs)
     if(box) graphics::box("plot")
     if(!missing(endexpr)) eval(parse(text=endexpr))
     if(write_t) title(sub=paste("t =", sprintf("%6.2f", t)), adj=1)
     if(sleep!=0) Sys.sleep(sleep)
     } # End for loop
  if(lastplot)
    {
    if(!missing(parexpr)) eval(parse(text=parexpr))
    plot(x, y, las=las, yaxt="n", type="n", log="y", xlim=xlim, ylim=ylim, ...)
    abline(h=lv$all, col=8)
    axis(2, lv$vals, lv$labs, las=las)
    do.call(points, args=owa(d=list(x=x, y=y), a=pointsarg, "x", "y"))
    if(box) graphics::box("plot")
    if(!missing(endexpr)) eval(parse(text=endexpr))
    }
  }
else if(log=="xy" | log=="yx") # -----------------------------------------------
  {
  lvx <- logVals(x, base=base)
  lvy <- logVals(y, base=base)
  for(t in allt)
     {
     if(!missing(parexpr)) eval(parse(text=parexpr))
     plot(x^(1/t), y^(1/t), las=las, xaxt="n", yaxt="n", type="n", xlim=xlim^(1/t), ylim=ylim^(1/t), ...)
     abline(h=(lvy$all)^(1/t), v=(lvx$all)^(1/t), col=8)
     axis(1, (lvx$vals)^(1/t), lvx$labs, las=las)
     axis(2, (lvy$vals)^(1/t), lvy$labs, las=las)
     pargs <- owa(d=list(x=x^(1/t), y=y^(1/t)), a=pointsarg, "x", "y")
     do.call(points, args=pargs)
     if(box) graphics::box("plot")
     if(!missing(endexpr)) eval(parse(text=endexpr))
     if(write_t) title(sub=paste("t =", sprintf("%6.2f", t)), adj=1)
     if(sleep!=0) Sys.sleep(sleep)
     } # End for loop
  if(lastplot)
    {
    if(!missing(parexpr)) eval(parse(text=parexpr))
    plot(x, y, las=las, xaxt="n", yaxt="n", type="n", log="xy", xlim=xlim, ylim=ylim, ...)
    abline(h=lvy$all, v=lvx$all, col=8)
    axis(1, lvx$vals, lvx$labs, las=las)
    axis(2, lvy$vals, lvy$labs, las=las)
    do.call(points, args=owa(d=list(x=x, y=y), a=pointsarg, "x", "y"))
    if(box) graphics::box("plot")
    if(!missing(endexpr)) eval(parse(text=endexpr))
    }
  }
else stop("log can only be 'x', 'y', or 'xy'.")
return(invisible(allt))
} # end of function ------------------------------------------------------------

# old way to get the transformation values:
#  if(length(steps)>1) steps <- max(steps, na.rm=TRUE)
#  # t-values for x^(1/t):
#  allt <- 10^(seq(0,2.5,len=1e4) )
#  # selection at upper half of these values;
#  # Otherwise, the animation slows down too much at the end
#  f <- 1.4 # multiplication factor due to length loss by unique
#  sel <- round(seq(1, 10, len=f*steps)^4)   #0.5*seq(1, 100, len=1.3*steps)^2 + 0.5*
#  sel2 <- unique(round(log10(seq(1, 10, len=f*steps))*f*steps))
#  sel2[1] <- 1
#  sel <- sel[sel2]
#  # final t-values for transition:
#  allt <- unique(round(allt[sel], 2))

# Current t value calculation is based on a polynomial of degree 17 fitted to
# the t values for 300 steps from this version

# very old way to get the transformation values:
# make more steps, as ca 35% are removed from allt later:
# steps <- round(steps * 1.35)#(1+(0.5+0.2*0.5)))
# steps <- 150
# allt <- 10^(seq(0,2.5,len=steps) )
# sel <- round(10^(seq(log10(steps/2), log10(steps), len=0.2*steps) ))
# allt <- allt[c(1:(sel[1]-1), sel)]

