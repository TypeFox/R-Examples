
#' @import animation

# TODO:
# density, stars, polygons? 
# plot3d - is this possible? persp is done...
# how to save while respecting intervals: interpolate data for frames
# in .do.loop:  
# EXAMPLE:
# ht <- 40
# x <- rep(1, ht)
# y <- 1:ht
# intervals <- jitter(ht:1/20) # in seconds
# times <- cumsum(intervals)
# ani.record(reset=T)
# 
# framerate <- 50 # per second
# intervals <- round(intervals*framerate)
# 
# # interpolate data? But, can I do this in general?
# # if so, the whole .do.loop would have to be rewritten to, indeed, smooth stuff.
# # a set of "smoothable" parameters.
# # if it were told to smooth, the intervals would be replaced by 1 and smoothed
# # versions of everything would be created. The resulting plot would be recorded.
# 
# # e.g. if intervals are 5,4,3,2,1. in [ x[1] , x[2] ) we want 5 intervals, etc.
# # in general for every variable: take x as sequence along it. it as y. then create xout
# spaced <- seq(min(times), max(times), length.out=max(times))
# x <- approx(times, x, xout=spaced)$y
# y <- approx(times, y, xout=spaced)$y

.setup.anim <- function (reset=TRUE, dev.control.enable=TRUE) {
  if (dev.cur()==1) dev.new()
  if (dev.control.enable) dev.control('enable')
  ani.record(reset=reset)
  # if (! is.null(interval)) .old.ani.options <<- ani.options(interval=interval)
}


.do.loop <- function(fn, times, show=TRUE, speed=1, use.times=TRUE, window=t,
      window.process=NULL, slice.args=list(), chunk.args=list(), oth.args=list(), 
      arg.dims=list(), chunkargs.ref.length=NULL) {
  # slice.args we take a slice and drop a dimension
  # chunk.args we cut without dropping
  # oth.args we leave alone
  mydiml <- function(obj) {
    if (is.null(dim(obj))) {
      if (length(obj)==1 || is.null(obj)) 0 else 1
    } else {
      length(dim(obj))
    }
  }
  
  utimes <- unique(times)
  nframes <- length(utimes)
  
  for (ar in names(slice.args)) {
    if (! ar %in% names(arg.dims)) arg.dims[[ar]] <- 0
    if (arg.dims[[ar]]==0) suppressWarnings(slice.args[[ar]] <- 
          rep(slice.args[[ar]], length=nframes))
  }
  if (! is.null(chunkargs.ref.length)) for (ar in names(chunk.args)) 
        suppressWarnings(chunk.args[[ar]] <- rep(chunk.args[[ar]], 
        length=chunkargs.ref.length))

  for (ca in names(chunk.args)) {
    chunk.args[[ca]] <- chunk.args[[ca]][order(times)]
  }
  times <- sort(times)
  utimes <- unique(times) # redo to make it correctly ordered
  
  mycalls <- list()
  for (t in 1:nframes) {
    # hack for anim.lines.formula
    win.t <- if (is.character(window)) eval(parse(text=window)) else eval(window)
    win.t <- win.t[win.t %in% 1:nframes]
    args.t <- list()
    for (an in names(slice.args)) {
      aa <- slice.args[[an]]
      dl <- mydiml(aa)
      args.t[[an]] <- if (dl <= arg.dims[[an]]) aa else switch(dl+1, aa, aa[t], 
            aa[,t], aa[,,t])
    }
    idx <- times %in% utimes[win.t]
    for (cn in names(chunk.args)) {
      ca <- chunk.args[[cn]]
      dl <- mydiml(ca)
      args.t[[cn]] <- switch(dl+1, ca, ca[idx], ca[,idx, drop=FALSE])
    }

    if (! is.null(window.process)) args.t <- window.process(args.t, times[idx])
    cl <- as.call(c(fn, args.t, oth.args)) # or match.call?
    mycalls[[t]] <- cl
  } 
  class(mycalls) <- "anim.frames"
  attr(mycalls, "speed") <- speed
  attr(mycalls, "times") <- if (use.times) utimes else order(utimes)
  attr(mycalls, "dev.control.enable") <- ! any(sapply(list(points, lines, text,
        symbols, segments, arrows), identical, fn))
  if (show) replay(mycalls)
  
  return(invisible(mycalls))
}

.col.interp <- function(colmat, smooth) {
  ncolours <- (ncol(colmat)-1)*smooth + 1
  colmat <- t(apply(colmat, 1, function(cl) 
    colorRampPalette(cl, alpha=TRUE)(ncolours)
  ))
  return(colmat) 
}

.interp <- function (obj, smooth) {
  size <- if(is.matrix(obj)) ncol(obj) else length(obj)
  xout <- seq(1, size, 1/smooth) 
  if (is.matrix(obj)) return(t(apply(obj, 1, function (y)
    approx(1:size, y, xout)$y)
  ))
  approx(1:size, obj, xout)$y
}

.plot.segments <- function(..., fn=quote(segments)) {
  mc <- match.call()
  dots <- list(...)
  if (! "xlab" %in% names(dots)) dots$xlab <- ""
  if (! "ylab" %in% names(dots)) dots$ylab <- ""
  plot(0,0,xlim=dots$xlim, ylim=dots$ylim, type="n", xlab=dots$xlab,
    ylab=dots$ylab, main=dots$main, sub=dots$sub)
  mc[[1]] <- fn
  mc$fn <- NULL
  eval(mc)
}

.plot.arrows <- function(...) .plot.segments(..., fn=quote(arrows))

#' Replay an \code{anim.frames} object
#' 
#' Replay all or some of the frames of an object.
#' 
#' @param x an \code{anim.frames} object
#' @param speed a new speed
#' @param frames numeric vector specifying which frames to replay
#' @param before an expression to evaluate before each frame is plotted
#' @param after an expression to evaluate after each frame is plotted
#' @param ... other arguments (not currently used)
#' 
#' @details
#' \code{before} and \code{after} will have the arguments from the
#' frame's call available in their environment - see the example.
#' 
#' The \code{plot} method simply calls \code{replay}.
#' 
#' @examples
#' 
#' myplot <- anim.plot(1:10, 1:10, speed=3)
#' replay(myplot, speed=5)
#' replay(myplot, frames=c(1,5,6,10))
#' 
#' myplot <- anim.plot(x<-rnorm(100), x+rnorm(100,0,3), 20, window=1:t, 
#'      show=FALSE, main="Regressions as sample size increases")
#' replay(myplot, after=abline(lm(y~x), col="red"))
#'  
#' @export
replay <- function(...) UseMethod("replay")

#' @export
#' @rdname replay
replay.anim.frames <- function(x, frames=1:length(x), speed=attr(x, "speed"),
  after=NULL, before=NULL, ...) {
  before2 <- substitute(before)
  after2 <- substitute(after)
  .setup.anim(dev.control.enable=attr(x, "dev.control.enable"))
  times <- attr(x, "times")
  intervals <- c(diff(times), 0)
  for (t in frames) {
    argl <- as.list(x[[t]])
    if (! missing(before)) eval(before2, argl)
    eval(x[[t]])
    if (! missing(after)) eval(after2, argl)
    ani.record()
    ani.pause(intervals[t]/speed)
  }
  invisible()
} 

#' @export
#' @rdname replay
plot.anim.frames <- function(x, ...) replay(x, ...)
  
#' Create an animated barplot.
#' 
#' @param height a vector, matrix or array. If a vector it is divided up by 
#'   \code{times} and \code{\link{barplot}} is called on each chunk. If a
#'   matrix, \code{\link{barplot}} is called on each column. If an array, 
#'   \code{\link{barplot}} is called on each matrix of form \code{height[,,i]}.
#' @param times a vector of times. If NULL and \code{height} is a matrix,
#'   the last dimension of \code{height} will be used
#' @param show,speed,use.times,window,window.process see \code{\link{anim.plot}} 
#' @param width,space,beside,names.arg,density,angle,col,border,horiz,xlim,ylim,xlab,ylab,main,sub,offset,legend.text,... arguments passed to \code{\link{barplot}}. 
#'  
#' @details
#' 
#' Arguments \code{width, names.arg, density, angle, col, border} 
#' and \code{offset} may be either vectors
#' of length \code{length(tbl)} or matrices with one column for each unique 
#' value of \code{times}. Other arguments should be length 1 or vectors.
#' 
#' @examples
#' anim.barplot(1:100, times=rep(1:10, each=10), ylim=c(0,100))
#' ## barplot with a matrix
#' ChickWeight$wq <- cut(ChickWeight$weight, 5)
#' tbl <- as.array(xtabs(~ wq + Diet + Time, data=ChickWeight))
#' ptbl <- prop.table(tbl, 2:3)
#' anim.barplot(ptbl, xlab="Diet", ylab="N", xlim=c(0,8), legend.text=paste(
#'      "Quintile", 1:5), col=1:5)
#' anim.barplot(tbl, xlab="Diet", ylab="N", beside=TRUE, ylim=c(0,20),
#'    legend.text=paste("Quintile", 1:5), col=1:5)
#'    
#' @export
anim.barplot <- function(...) UseMethod("anim.barplot")

#' @export
#' @rdname anim.barplot
anim.barplot.default <- function(height, times=NULL, 
      show=TRUE, speed=1, use.times=TRUE, window=t, window.process=NULL, 
      width=1, space=NULL, names.arg=NULL, beside=FALSE, density=NULL, 
      angle=NULL, col=NULL, border=NULL, horiz=FALSE, xlim=NULL, 
      ylim=NULL, xlab=NULL, ylab=NULL, main=NULL, sub=NULL, offset=NULL, 
      legend.text=NULL, ...) {
  # plot data
  slice.args <- list(height=height, space=space, xlim=xlim, ylim=ylim, main=main, 
        sub=sub, xlab=xlab, ylab=ylab, legend.text=legend.text, width=width, 
        names.arg=names.arg, density=density, angle=angle, border=border, 
        offset=offset, col=col)

  args <- list(...)  
  ltdim <- if (is.logical(legend.text)) 0 else 1
  oth.args <- args
  oth.args$beside <- beside
  chunk.args <- list()
  if (is.vector(height)) chunk.args$height=height else slice.args$height=height

  hdim <- if(is.matrix(height)) 1 else 2
  if (is.null(times)) {
    if (is.array(height)) times <- 1:tail(dim(height), 1) else 
          stop("'times' not specified")
  } else if (length(times)==1) {
    lng <- if (is.array(height)) tail(dim(height), 1) else length(height)
    if (lng %% times != 0) warning("'height' length is not an exact multiple of 'times'")
    times <- rep(1:times, each=lng/times)
  }
  crl <- if(is.vector(height)) max(length(height), length(times))

  arg.dims <- list(height=hdim, space=1, xlim=1, ylim=1, main=0, sub=0, xlab=0, 
        ylab=0, space=1, legend.text=ltdim, col=1, density=1, angle=1, names.arg=1,
        border=1, offset=1, width=1)
  .do.loop(barplot, times=times, use.times=use.times, window=substitute(window),
        window.process=window.process, show=show, speed=speed, 
        slice.args=slice.args, chunk.args=chunk.args, 
        oth.args=oth.args, arg.dims=arg.dims, chunkargs.ref.length=crl)
}

# gapminder data: life expectancy, regions (colour), GDP/cap ppp infadjust, pop total,
# years

#' Create an animated plot.
#' 
#' \code{anim.plot}
#' 
#' @param x,y vectors of x and y coordinates. These can be passed in any way 
#'   accepted by \code{\link{xy.coords}}.
#' @param times a vector of times. If \code{times} is length one, there will
#'   be that many frames, equally divided over the length of \code{x} and \code{y}.
#' @param show if false, do not show plot; just return calls.
#' @param speed animation speed.
#' @param window what window of times to show in each animation. The default,
#'   \code{t}, shows just plots from time t. To draw a plot incrementally,
#'   use \code{window=1:t}. 
#' @param window.process function to call on each window of each times. See details.
#' @param use.times if \code{TRUE}, animation speed is determined by the 
#'   \code{times} argument. If \code{FALSE}, animation speed is constant.
#' @param xlim,ylim,col,pch,labels,cex,lty,lwd,asp,xaxp,yaxp,... arguments passed to 
#'   \code{\link{plot}}.
#' @param fn function called to create each frame
#' @param data a data frame from where the values in \code{formula} should be 
#'    taken
#' @param formula a \code{\link{formula}} such as \code{y ~ x + time}
#' @param subset a vector specifying which rows of \code{data} to use
#'   
#' @details
#' 
#' Each unique level of \code{times} will generate a single frame of animation. 
#' The frames will be ordered by \code{times}.
#' 
#' In general:
#' 
#' \itemize{ 
#' \item Parameters that apply to each point of the plot, such as 
#' \code{xlim, ylim, col, pch, labels} and \code{cex}, can be passed as vectors 
#' which will be recycled to \code{length(times)}. 
#' \item Parameters that apply
#' to the plot as a whole, and always have length 1, such as \code{xlab} and
#' \code{main}, can be passed as vectors and will be recycled to the number of
#' frames. 
#' \item Parameters that apply to the plot as a whole, and can have
#' length > 1, such as \code{xlim} and \code{ylim}, can be passed as vectors or
#' matrices. If vectors, the same vector will be passed to every frame. If
#' matrices, column \code{i} will be passed to the \code{i}'th frame. 
#' }
#' 
#' \code{window.process} should be a function which takes
#' two arguments: a list of potential arguments for the underlying
#' call to \code{plot}, and a vector of times. The function should return
#' the list of arguments after modification. This allows e.g. drawing 
#' "trails" of plot points. See the example
#' 
#' @examples
#' x <- rep(1:100/10, 10)
#' times <- rep(1:10, each=100)
#' y <- sin(x*times/4)
#' anim.plot(x,y,times, type="l", col="orange", lwd=2)
#' 
#' ## changing colours - a per-point parameter
#' anim.plot(x,y,times, ylab="Sine wave", type="p", col=rainbow(100)[x *10])
#' 
#' ## changing line width - a whole-plot parameter
#' anim.plot(x, y, times, lwd=1:10, type="l")
#'           
#' ## times as a single number 
#' anim.plot(1:10, 1:10, times=5)
#'            
#' ## incremental plot
#' anim.plot(1:10, 1:10, window=1:t)
#' 
#' ## moving window
#' anim.plot(1:10, 1:10, window=(t-2):t)
#' 
#' ## Formula interface
#' ChickWeight$chn <- as.numeric(as.factor(ChickWeight$Chick))
#' tmp <- anim.plot(weight ~ chn + Time, data=ChickWeight, col=as.numeric(Diet), 
#'      pch=as.numeric(Diet), speed=3)
#' 
#' # adding extra arguments:
#' replay(tmp, after=legend("topleft", legend=paste("Diet", 1:4), pch=1:4, col=1:4))
#'  
#'  ## Zooming in:
#'  x <- rnorm(4000); y<- rnorm(4000)
#'  x <- rep(x, 10); y <- rep(y, 10)
#'  xlims <- 4*2^(-(1:10/10))
#'  ylims <- xlims <- rbind(xlims, -xlims) 
#'  anim.plot(x, y, times=10, speed=5, xlim=xlims, ylim=ylims, 
#'        col=rgb(0,0,0,.3), pch=19)
#'  
#'  ## window.process to create a faded "trail":
#'  anim.plot(1:50, 1:50, speed=12, window=t:(t+5), 
#'        window.process=function(args, times){
#'          times <- times - min(times) 
#'          alpha <- times/max(times)
#'          alpha[is.na(alpha)] <- 1
#'          args$col <- rgb(0,0,0, alpha)
#'          return(args)
#'        })
#'        
#'  ## gapminder plot:
#'  pl <- palette(adjustcolor(rainbow(23), 1, .6, .6, .6, 
#'        offset=c(0,0,0,-0.1)))
#'  anim.plot(lifex ~ GDP + year, data=gm_data, log="x", 
#'       cex=sqrt(pop)*0.0004, pch=19, col=region, xlab="GDP", 
#'       ylab="Life expectancy", speed=10, subset=year > 1850 & !year %% 5)
#'  palette(pl)
#'  
#'  \dontrun{
#'  ## Earthquakes this week
#'  if (require('maps')) {
#'    eq = read.table(
#'        file="http://earthquake.usgs.gov/earthquakes/catalogs/eqs7day-M1.txt", 
#'        fill=TRUE, sep=",", header=TRUE)
#'    eq$time <- as.numeric(strptime(eq$Datetime, "%A, %B %d, %Y %X UTC"))
#'  eq <- eq[-1,]
#'    map('world')
#'    maxdepth <- max(max(eq$Depth), 200)
#'    tmp <- anim.points(Lat ~ Lon + time, data=eq, cex=Magnitude, col=rgb(
#'          1-Depth/maxdepth, 0, Depth/maxdepth,.7), pch=19, speed=3600*12, 
#'          show=FALSE)   
#'    replay(tmp, before=map('world', fill=TRUE, col="wheat"))
#'  }
#'  
#'  
#'  ## Minard's plot
#'  if (require('maps')) {
#'    map('world', xlim=c(22, 40), ylim=c(52,58))
#'    title("March of the Grande Armee on Moscow")
#'    points(cities$long, cities$lat, pch=18)
#'    text(cities$long, cities$lat, labels=cities$city, pos=4, cex=.7)
#'    with(troops[troops$group==1,], anim.lines(x=long, 
#'          y=lat, window=t:(t+1), speed=3, lwd=survivors/10000))
#' }
#' }
#' @export
anim.plot <- function(...) UseMethod("anim.plot")

#' @export
#' @rdname anim.plot
anim.points <- function(...) UseMethod("anim.points")

#' @export
#' @rdname anim.plot
anim.lines <-function(...) UseMethod("anim.lines")

#' @export
#' @rdname anim.plot
anim.text <-function(...) UseMethod("anim.text")

#' @export 
#' @rdname anim.plot
anim.plot.default <- function (x, y=NULL, times=1:length(x), speed=1, show=TRUE,
      use.times=TRUE, window=if (identical(fn, lines)) t:(t+1) else t, window.process=NULL, xlim=NULL, ylim=NULL, 
      col=par("col"), xaxp=NULL, yaxp=NULL, pch=par("pch"), cex=1, labels=NULL, 
      asp=NULL, lty=par("lty"), lwd=par("lwd"), fn=plot, ...) {  
  
  args <- list(...)
  if (! "xlab" %in% names(args)) args$xlab <- deparse(substitute(x))
  if (! "ylab" %in% names(args)) args$ylab <- deparse(substitute(y))
  xy <- xy.coords(x, y, recycle=TRUE)
  x <- xy$x
  y <- xy$y
  args$xlim <- if (is.null(xlim)) range(x[is.finite(x)]) else xlim
  args$ylim <- if (is.null(ylim)) range(y[is.finite(y)]) else ylim
  
  if (length(times)==1) {
    lng <- length(x)
    if (lng %% times != 0) warning("'height' length is not an exact multiple of 'times'")
    times <- rep(1:times, each=lng/times)
  }
  chunk.args <- list(x=x, y=y, col=col, pch=pch, cex=cex)
  slice.args <- c(list(asp=asp, lty=lty, lwd=lwd, xaxp=xaxp, yaxp=yaxp), args)
  if (identical(fn, text)) {
    chunk.args$labels <- labels
    slice.args$labels <- NULL
  }
  .do.loop(fn, times=times, speed=speed, show=show, use.times=use.times, 
        window=substitute(window), window.process=window.process, 
        chunk.args=chunk.args, slice.args=slice.args, 
        arg.dims=list(xlab=0, ylab=0, xlim=1, ylim=1, xaxp=1, yaxp=1, lwd=0, 
        lty=0, asp=0, panel.first=1, panel.last=1, x=1, y=1, col=1, pch=1, cex=1, 
        type=0), chunkargs.ref.length=max(length(x), length(y)))
}

#' @export 
#' @rdname anim.plot
anim.plot.formula <- function(formula, data=parent.frame(), subset=NULL, 
      fn=plot, window=t, ...) {
  if (missing(formula) || ! inherits(formula, "formula")) 
    stop("'formula' missing or invalid")
  
  # cargo-culted from plot.formula
  m <- match.call(expand.dots=FALSE)
  eframe <- parent.frame() 
  md <- eval(m$data, eframe)
  dots <- lapply(m$..., eval, md, eframe) 
  mf <- model.frame(formula, data=md)
  subset.expr <- m$subset
  if (!missing(subset)) {
    s <- eval(subset.expr, data, eframe)
    l <- nrow(mf)
    dosub <- function(x) if (length(x) == l) x[s] else x
    dots <- lapply(dots, dosub)
    mf <- mf[s, ]
  }
  
  # get levels of t. 
  x <- mf[,2]
  y <- mf[,1]
  tm <- if (ncol(mf) >= 3) mf[,3] else 1:length(x)
  
  # why doesn't ordering happen happen OK in .do.loop?
  ot <- order(tm)
  # we are basically praying here:
  dots <- lapply(dots, function(z) if (length(z)==length(tm)) z[ot] else z) 
  
  x <- x[ot]
  y <- y[ot]
  tm <- tm[ot]
  if (! "xlab" %in% names(dots)) dots$xlab <- all.vars(mf)[2] 
  if (! "ylab" %in% names(dots)) dots$ylab <- all.vars(mf)[1]
  do.call("anim.plot", c(list(x=x, y=y, times=tm, window=substitute(window), fn=fn), dots))
}

#' @export 
#' @rdname anim.plot
anim.points.default <- function(...) anim.plot.default(..., fn=points)

#' @export 
#' @rdname anim.plot
anim.lines.default <- function(...) anim.plot.default(..., fn=lines)

#' @export 
#' @rdname anim.plot
anim.text.default <- function(...) anim.plot.default(..., fn=text)

#' @export 
#' @rdname anim.plot
anim.symbols <- function(...) anim.plot.default(..., fn=symbols)


#' @export 
#' @rdname anim.plot
anim.points.formula <- function(formula, ...) {
  m <- match.call(expand.dots=TRUE)
  fn <- as.character(m[[1]])
  fn <- sub("anim.([a-z]+).formula", "\\1", fn)
  fn <- eval(as.name(fn))
  m[[1]] <- quote(anim.plot.formula)
  m[["fn"]] <- fn
  eval(m)
}

#' @export 
#' @rdname anim.plot
anim.lines.formula <- anim.points.formula

#' @export 
#' @rdname anim.plot
anim.text.formula <- anim.points.formula


#' Create an animated contour plot or perspective plot
#' 
#' Create an animated contour plot or perspective plot of 3D data.
#' 
#' @param x,y,z,... arguments passed to \code{\link{contour}} or \code{\link{persp}}
#' @param times,speed,use.times,window,window.process,show see 
#'    \code{\link{anim.plot}} for details.
#' @param fn underlying function to use.
#' 
#' @examples
#' 
#' tmp <- volcano
#' tmp[] <- 200 - ((row(tmp) - 43)^2 + (col(tmp) - 30)^2)/20
#' cplot <- array(NA, dim=c(87,61,20))
#' cplot[,,1] <- tmp
#' cplot[,,20] <- volcano
#' cplot <- apply(cplot, 1:2, function(x) seq(x[1], x[20], length.out=20))
#' cplot <- aperm(cplot, c(2,3,1))
#' anim.contour(z=cplot, times=1:20, speed=3, levels=80 + 1:12*10, lty=c(1,2,2))
#' anim.filled.contour(z=cplot, times=1:20, speed=3, levels=80 + 1:12*10, 
#'    color.palette=terrain.colors)
#'    
#' cplot2 <- apply(cplot, 1:2, function(x) seq(0, x[20], length.out=20))
#' cplot2 <- aperm(cplot2, c(2,3,1))
#' anim.persp(z=cplot2, times=1:20, xlab="", ylab="", zlab="Height", phi=45,
#' theta=30, speed=5, border=NA, r=3, col="yellowgreen", shade=.5, box=FALSE)
#'  
#' @export
anim.contour <- function(...) UseMethod("anim.contour")

#' @export
#' @rdname anim.contour
anim.filled.contour <- function(...) UseMethod("anim.filled.contour")

#' @export
#' @rdname anim.contour
anim.filled.contour.default <- function(...) anim.contour.default(..., fn=filled.contour)

#' @export
#' @rdname anim.contour
anim.persp <- function(...) {
  m <- match.call(expand.dots=TRUE)
  m[[1]] <- quote(anim.contour)
  m$fn <- quote(persp)
  eval(m)
}

#' @export
#' @rdname anim.contour
anim.contour.default <- function(x, y, z, times, speed=1, use.times=TRUE, window=t, 
      window.process=NULL, show=TRUE, fn=contour, ...) {
  if (missing(z)) {
    z <- x 
    x <- seq(0,1, length.out=dim(z)[1])
    y <- seq(0,1, length.out=dim(z)[2])
  }
  if (missing(x)) x <- seq(0,1, length.out=dim(z)[1])
  if (missing(y)) y <- seq(0,1, length.out=dim(z)[2])
  dots <- list(...)
  slice.args <- list(z=z)
  slice.args$x <- x
  slice.args$y <- y
  if (! "zlim" %in% names(dots)) dots$zlim <- range(z, finite=TRUE)
  if (! "xlim" %in% names(dots)) dots$xlim <- range(x, finite=TRUE)
  if (! "ylim" %in% names(dots)) dots$ylim <- range(y, finite=TRUE)
  if (length(times)==1) times <- 1:times
  .do.loop(fn, times=times, show=show, use.times=use.times,
        window=substitute(window), window.process=window.process,
        slice.args=c(slice.args, dots), 
        arg.dims=list(z=2, x=1, y=1, nlevels=0, levels=1, 
        labels=1, labcex=0, drawlabels=0, xlim=1, ylim=1, zlim=1, vfont=1,
        axes=0, frame.plot=0, col=1, lty=1, lwd=1, color.palette=1))
}

#' Draw an animated histogram.
#' 
#' @param x,density,angle,col,border,... parameters passed to \code{\link{hist}}.
#' @param times,show,speed,use.times,window,window.process see 
#'    \code{\link{anim.plot}}.
#' 
#' @details
#' Parameters \code{x, density, angle, col} and \code{border} are all
#' "chunked", i.e. first recycled to the length of \code{times} or \code{x}
#' (whichever is longer), then split according to the unique values of \code{times}.
#' See \code{\link{anim.plot}} for more details.
#' 
#' @examples 
#' anim.hist(rep(rnorm(5000), 7), times=rep(1:7, each=5000), 
#'      breaks=c(5,10,20,50,100,200, 500, 1000))
#' @export
anim.hist <- function(x, times, speed=1, show=TRUE, use.times=TRUE, window=t,
      window.process=NULL, density=NULL, angle=NULL, col=NULL, border=NULL, ...) {
  
  dots <- list(...)
  if (! "breaks" %in% names(dots)) dots$breaks = "Sturges"
  if (! "xlab" %in% names(dots)) dots$xlab <- ""
  if (! "main" %in% names(dots)) dots$main <- "Histogram"
    
  dbr <- if (is.matrix(dots$breaks)) 1 else 0
  .do.loop(hist, times=times, show=show, speed=speed, use.times=use.times, 
        window=substitute(window), window.process=window.process, 
        chunk.args=list(x=x, density=density, angle=angle, col=col, 
        border=border), slice.args=dots, arg.dims=list(breaks=dbr, xlim=1, 
        ylim=1, xlab=1, x=1), chunkargs.ref.length=max(length(x), length(times)))
}


#' Draw an animation of line segments or arrows.
#' 
#' @param x0,y0,x1,y1,col,lty,lwd,length,angle,code,... arguments passed to \code{\link{segments}} or
#'     \code{\link{arrows}}
#' @param times,speed,show,use.times,window,window.process see \code{\link{anim.plot}} for details
#' @param fn underlying function to use
#' 
#' @details
#' 
#' \code{anim.segments} and \code{anim.arrows} draw lines on to an existing plot.
#' If you want to redraw the plot between each frame, use \code{anim.arrowplot}
#' or \code{anim.segmentplot}.
#' 
#' If both \code{x1} and \code{y1} are missing, then segments are plotted
#' from the current time to the following time in each frame. If only \code{x1}
#' is missing it is set equal to \code{x0}, similarly if only \code{y1} is 
#' missing.
#' 
#' @examples
#' anim.segments(x0=rep(1:5, 5), y0=rep(1:5, each=5), y1=rep(2:6, each=5), 
#'      times=rep(1:5, each=5) )
#'  
#' ## Short version
#' anim.arrowplot(rep(1:5, 5), rep(1:5, each=5), times=5)
#' 
#' if (require('maps')) {
#'    hr <- subset(hurricanes, lat > 0 & lat < 50 & lon > -95 & lon < -20 & 
#'          Shour %% 6 == 0)
#'    hr$dlat <- cos(hr$diruv/360*2*pi) * hr$maguv / 8
#'    hr$dlon <- sin(hr$diruv/360*2*pi) * hr$maguv / 8
#'    hr$name <- sub("\\s+$", "", hr$name)
#'    map('world', xlim=c(-95,-20), ylim=c(0,50))
#'    title("Hurricanes, 2009")
#'    with(hr[!duplicated(hr$name),], text(lon, lat, 
#'          labels=paste0(name, "\n", Yr), cex=0.8))
#'    with(hr, anim.arrows(x0=lon, y0=lat, y1=lat+dlat, x1=lon+dlon, 
#'          times=Shour, speed=12, col=rgb(0,0,1,0.8), length=.1, lwd=2)) 
#' }
#' @export
anim.segments <- function(x0, y0, x1=NULL, y1=NULL, times=NULL, speed=1, show=TRUE, 
      use.times=TRUE, window=t, window.process=NULL, fn=segments, 
      col=NULL, lty=NULL, lwd=NULL, ...) {
  dots <- list(...)
  if (! "xlim" %in% names(dots)) dots$xlim <- range(c(x0, x1), na.rm=T)
  if (! "ylim" %in% names(dots)) dots$ylim <- range(c(y0, y1), na.rm=T)
  
  crl <- max(length(x0), length(x1), length(y0), length(y1), na.rm=T)
  if (is.null(times)) times <- 1:crl
  if (length(times)==1) {
    if (crl %% times != 0) warning(
      "length of longest vector is not an exact multiple of 'times'")
    times <- rep(1:times, each=crl/times)
  }
  
  if (is.null(x1) && is.null(y1)) {
    x1 <- x0[times > min(times)]
    x0 <- x0[times < max(times)]
    y1 <- y0[times > min(times)]
    y0 <- y0[times < max(times)]
    times <- times[times > min(times)]
  } else if (is.null(x1)) x1 <- x0 else if (is.null(y1)) y1 <- y0
  
  chunk.args <- list(x0=x0, y0=y0, x1=x1, y1=y1, col=col, lty=lty, lwd=lwd)
  for (ca in c("length", "angle", "code")) if (ca %in% names(dots)) 
        chunk.args[[ca]] <- dots[[ca]]
  .do.loop(fn, times=times, show=show, speed=speed, use.times=use.times, 
        window=substitute(window), window.process=window.process, 
        chunk.args=chunk.args,   
        slice.args=dots, arg.dims=list(xlim=1, ylim=1), chunkargs.ref.length=crl)
}


#' @export
#' @rdname anim.segments
anim.arrows <- function(..., length=0.25, angle=30, code=2) anim.segments(...,
      length=length, angle=angle, code=code, fn=arrows)

#' @export
#' @rdname anim.segments
anim.segmentplot <- function(...) anim.segments(..., 
      fn=.plot.segments)


#' @export
#' @rdname anim.segments
anim.arrowplot <- function(...) anim.segments(..., 
      fn=.plot.arrows)

#' Draw an animated curve.
#' 
#' This function is the animated version of \code{\link{curve}}.
#' 
#' @param expr a function which takes two arguments, or an expression involving
#'    \code{x} and \code{t}.
#' @param x values of \code{x} at which the function will be evaluated in each frame.
#'    Alternatively, you may specify \code{from, to} and \code{n}.
#' @param from,to endpoints of \code{x}
#' @param n number of values of \code{x} at which the function will be evaluated
#'   for each frame
#' @param times vector of values of \code{t} at which the function will be 
#'   evaluated. Each unique value creates a single animation frame.
#' @param type,... parameters passed to \code{\link{anim.plot.default}}
#' 
#' @details
#' Note that \code{times} is interpreted differently here than elsewhere. In
#' particular, it cannot be a length-1 vector.
#' 
#' @examples
#' anim.curve(x^t, times=10:50/10, n=20)
#' anim.curve(sin(x*t), times=1:30, n=100, speed=12, col="darkgreen", from=-1, to=1)
#' 
#' ## curve is constant in t, but window moves. 
#' ## NB: 'from' and 'to' control where the expression is evaluated. 
#' ## 'xlim' just controls the window.
#' anim.curve(sin(cos(-x)*exp(x/2)), times=0:100/10, from=-5, to=10, n=500, 
#'      col="red", lwd=2, xlim=rbind(top <- seq(-5, 10, 1/10), top+5))
#' @export
anim.curve <- function(expr, x=NULL, from=0, to=1, n=255, times, type="l", ...) {
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    expr <- call(as.character(sexpr), as.name("x"), as.name("t"))
  } else {
    expr <- sexpr
  }
  if (is.null(x)) x <- seq.int(from, to, length.out=n)
  
  y <- outer(x, times, function (x,t) {
    ll <- list(x=x, t=t)
    eval(expr, envir=ll, enclos=parent.frame())
  })
  y <- as.vector(y)
  times <- rep(times, each=length(x))
  anim.plot(x=x, y=y, times=times, type=type, ...)
}

#' Save an anim.frames object in various formats.
#' 
#' This function simply calls replay on the object and then calls
#' \code{\link{saveGIF}} and friends on the result.
#' 
#' @param obj an \code{anim.frames} object
#' @param type one of 'GIF', 'Video', 'SWF', 'HTML', or 'Latex'
#' @param filename file to save to
#' @param ... arguments passed to e.g. \code{\link{saveGIF}}
#' 
#' @details
#' 
#' For most of the underlying functions, the \code{times} parameter
#' will be ignored as if you had set \code{use.times=FALSE}.
#' 
#' @examples
#' 
#' \dontrun{
#' tmp <- anim.plot(1:10, 1:10, pch=1:10, show=FALSE)
#' anim.save(tmp, "GIF", "filename.gif")
#' 
#' ## for anything more complex. Note the curlies:
#' saveGIF({replay(tmp, after=legend("topleft", legend="My legend"))},
#'  "filename.gif")
#' }
#' @export
anim.save <- function(obj, type, filename, ...) {
  stopifnot(type %in% c("GIF", "Video", "SWF", "HTML", "Latex"))
  fn <- as.name(paste("save", type, sep=""))
  mf <- match.call(expand.dots=FALSE)
  mf[[1]] <- fn
  mf$obj <- NULL
  mf$expr <- substitute(replay(obj))
  mf$type <- NULL
  switch(type, 
    "GIF"=mf$movie.name <- filename, 
    "Video"=mf$video.name <- filename,
    "SWF"=mf$swf.name <- filename, 
    "HTML"=mf$htmlfile <- filename, 
    "Latex"=mf$latex.filename <- filename
  )
  eval(mf)
}

#' Merge anim.frames objects
#' 
#' Merge two or more anim.frames objects to create a new anim.frames object
#' 
#' @param ... anim.frames objects returned from, e.g. \code{\link{anim.plot}}
#' @param speed speed for the merged object. This may be left unspecified only
#'    if all objects have the same speed.
#'
#' @details
#' If two or more calls in the merged animation are at the same time, calls
#' from the earlier object in \code{...} will be run first. 
#' 
#' If you merge two animations from \code{\link{anim.plot}}, plot.window will be
#' called before each frame of the merged animation. This may not be what
#' you want. Instead, use \code{anim.points} or similar for all but the first
#' animation.
#' 
#' 
#' @examples
#' tmp <- anim.plot(1:5, 1:5, speed=2)
#' tmp2 <- anim.plot(1:5, 5:1, col="red", speed=2)
#' ## Not what you want:
#' replay(merge(tmp, tmp2))
#' 
#' ## better:
#' tmp3 <- anim.points(1:5, 5:1,col="red", speed=2)
#' newf <- merge(tmp, tmp3)
#' replay(newf)
#' ## NB: result of the merge looks different from the two
#' ## individual animations
#' 
#' ## not the same:
#' newf2 <- merge(tmp2, tmp) 
#' ## points will be called before plot!
#' replay(newf2)
#' @export
merge.anim.frames <- function(..., speed=NULL) {
  frs <- list(...)
  speeds <- sapply(frs, attr, "speed")
  if (is.null(speed)) if (max(speeds) > min(speeds)) {
    stop("'speed' not specified but some objects have different speeds")
  } else {
    speed <- max(speeds)
  }
  times <- c(sapply(frs, attr, "times"))
  tiebreaks <- sapply(1:length(frs), function(x) rep(x, length(frs[[x]])))
  newfr <- unlist(frs, recursive=FALSE)
  newfr <- newfr[order(times, tiebreaks)]
  times <- sort(times)
  class(newfr) <- "anim.frames"
  attr(newfr, "times") <- times
  attr(newfr, "speed") <- speed
  # maybe this is the wrong way to think about it?
  attr(newfr, "dev.control.enable") <- any(sapply(frs, attr, "dev.control.enable"))
  newfr
}



#' Troop numbers for the Grande Armee's march on Moscow
#' @name troops
NULL

#' Cities near the Grande Armee's march on Moscow
#' @name cities
NULL

#' Temperatures for the Grande Armee's march on Moscow
#' @name temps
NULL


#' Wind speed data for hurricanes in 2009
#' @name hurricanes
#' @source http://myweb.fsu.edu/jelsner/Data.html
NULL

#' Gapminder GDP, life expectancy and population data
#' @name gm_data
#' @source http://gapminder.org
NULL
