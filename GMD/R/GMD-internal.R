
## ------------------------------------------------------------------------
## Wrapping the Library Routine
## ------------------------------------------------------------------------
.mdpa <- function(v1, v2){
  if (length(v1) != length(v2)){
    stop("The `mdpa` function should take two equal-length vectors as input!\n")
  }
  if (sum(v1) != sum(v2)){
    stop("The `mdpa` function should take two equal-mass vectors as input!\n")
  }
  .Call("mdpa", v1, v2, PACKAGE="GMD")
}


.gmd0 <- function(v1, v2, pseudocount=0){
  if (length(v1) != length(v2)){
    stop("The `gmd0` function should take two equal-length vectors as input!\n")
  }
  res <- .Call("gmd0", v1, v2, pseudocount, PACKAGE="GMD")
  return(res)
}


.gmd <- function(v1, v2, pseudocount=0){
  res <- .Call("gmd", v1, v2, pseudocount, PACKAGE="GMD")
  return(res)
}


## ------------------------------------------------------------------------
## Other internal functions
## ------------------------------------------------------------------------
.wordwrap <-
  function(x,len)
{
  l <- nchar(x)
  m <- matrix(ncol=2,nrow=ceiling(l/len))
  for (i in 1:nrow(m)){
    m[i,] <- c(1+len*(i-1),min(len*i,l))
  }
  res <- mapply(FUN=substr,m[,1],m[,2],MoreArgs=list(x=x))
  paste(res,sep="",collapse="\n")
}


.is.gmd <-
  function(x)
{
  ## ##"gmd" %in% class(x)
  inherits(x,"gmd")
}


.is.gmdp <-
  function(x)
{
  ## ##"gmdp" %in% class(x)
  inherits(x,"gmdp")
}


.is.gmdm <-
  function(x)
{
  ## ##"gmdm" %in% class(x)
  inherits(x,"gmdm")
}


.is.grouped <-
  function(x)
{
  x <- as.character(x)
  x.f <- factor(x,levels=unique(x),ordered=TRUE)
  identical(as.character(sort(x.f)),x)
}


.resolveBin <-
  function(l,r)
{
  lower.bound <- (1:ceiling(l/r)-1)*r+1
  upper.bound <- (1:ceiling(l/r))*r
  upper.bound <- sapply(upper.bound, function(x) ifelse(x>l,l,x))
  cbind(lower.bound,upper.bound)
}


.resolveHist <-
  function(v,r)
{
  l <- length(v)
  b <- .resolveBin(l,r)
  sapply(1:nrow(b),function(i) sum(v[b[i,1]:b[i,2]]))
}


##colorpanel(16,"gray96","gray47")
.grays16 <-
    c(
      "#F5F5F5", "#EDEDED", "#E4E4E4", "#DCDCDC", "#D4D4D4", "#CBCBCB", "#C3C3C3", "#BBBBBB",
      "#B2B2B2", "#AAAAAA", "#A2A2A2", "#999999", "#919191", "#898989", "#808080", "#787878"
      )


.setTextContrastColor <- function(color){
  ifelse( mean(col2rgb(color)) > 127, "black", "white")
}


##' A set of colors generated from palette Dark2
##'
##' This is a realization of `brewer.pal(n, name="Dark2")' function in the Package `RColorBrewer'.
##' @title A set of colors generated from palette Dark2
##' @param n integer, indicating the number of colors to generate
##' @return a set of colors
.brewer.pal.Dark2 <- function(n){
  tmp <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  if(n>8){
    warning("`.brewer.pal.Dark2' can generate no more than 8 unique colors; colors are recycled.")
    colors <- rep(tmp,ceiling(n/length(tmp)))[1:n]
  } else {
    colors <- tmp[1:n]
  }
  colors
}


.scale.range <-
  function(r,f)
{
  range(c(r/f,r*f))
}


##' Scale values to a new range: c(low, high)
##'
##' Scale values to a new range: c(low, high)
##' @title Scale values to a new range.
##' @param x numeric
##' @param low numeric, lower bound of target values
##' @param high numeric, higher bound of target values
##' @return an object with the same dimention of `x'.
.scale.x <-
  function(x,low=0,high=1,na.rm=TRUE)
{
  if(identical(max(x,na.rm=na.rm),min(x,na.rm=na.rm))) NA
  a <- 1/(max(x)-min(x))
  b <- -min(x)/(max(x)-min(x))
  a*x+b
}



##' Scale values to make them follow Standard Normal Distribution
##'
##' Scale values to make them follow Standard Normal Distribution
##' @title Scale values to make them follow Standard Normal Distribution
##' @param x numeric
##' @param scale character, indicating the type to scale.
##' @param na.rm logical
##' @return an object with the same dimention of `x'.
.scale.data <-
  function(x,scale,na.rm=TRUE)
{
  if(scale=="row"){
    x <- sweep(x,1,rowMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,1,sd,na.rm=na.rm)
    x <- sweep(x,1,sx,FUN="/")
  } else if(scale=="column"){
    x <- sweep(x,2,colMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,2,sd,na.rm=na.rm)
    x <- sweep(x,2,sx,,FUN="/")
  }
  x
}


##' Call a function with arguments
##'
##' Call a function with arguments
##' @title Call a function with arguments
##' @param FUN function or function name
##' @param ... unnameed function arguments
##' @param MoreArgs named (or unnameed) function arguments
.call.FUN <-
  function(FUN,...,MoreArgs)
{
  FUN <- match.fun(FUN)
  tmp.MoreArgs <- list(...)
  if (!.invalid(MoreArgs)){
    if (length(MoreArgs)>=1){
      for (i in 1:length(MoreArgs)) tmp.MoreArgs[[names(MoreArgs)[i]]] <- MoreArgs[[i]]
    }
  }
  ret <- do.call(FUN, tmp.MoreArgs)
  
##   attr(ret,"call") <-
##     sprintf("%s%s",
##             as.character(quote(FUN)),
##             substring(capture.output(dput(MoreArgs,control=c("keepNA", "keepInteger"))),5)
##             )
  
  if ("call" %in% names(ret)){
    ret$call <- match.call()
  }
  if ("call" %in% names(attributes(ret))){
    attr(ret,"call") <- match.call()
  }
  return(ret)
}


##' Plot text 
##'
##' Plot text 
##' @title Plot text 
##' @param x character, text to plot
##' @param cex 
##' @param forecolor color of foreground
##' @param bg color of background
##' @param bordercolor color of border
##' @param axes as in \code{graphics:::plot}
##' @param ... additional arguments for \code{graphics:::text}
.plot.text <- function(x,xlim=c(0,1),ylim=c(0,1),cex=1,forecolor=par("fg"),bg=par("bg"),bordercolor=NA,axes=FALSE,...){
  if (.invalid(x)){
    x <- NULL
  }
  if (is.null(x)){
    x <- ""
  } else if (is.na(x)){
    x <- 'NA'
  }
  
  plot(xlim,ylim,type="n",ylab="",xlab="",xaxt="n",yaxt="n",axes=axes)
  rect(xleft=0, ybottom=0, xright=1, ytop=1, col=bg, border=bordercolor)
  text(0.5,0.5,x,cex=cex,...)
}




##' Normalize/scale a vector given both source and target ranges
##'
##' Normalize/scale a vector given both source and target ranges
##' @title Normalize/scale a vector given both source and target ranges
##' @param v a numeric vector.
##' @param source.x1 the lower bound of the original vector; the default is the minimal value of \code{v}.
##' @param source.x2 the upper bound of the original vector; the default is the maximal value of \code{v}.
##' @param target.x1 the lower bound of the target vector; the default is \code{0}.
##' @param target.x2 the upper bound of the target vector; the default is \code{1}.
##' @param na.rm a logical indicating whether missing values should be removed; \code{TRUE} by default.
##' @return a vector 
##' @examples
##' v <- 1:10
##' normalizeVector(v,NULL,NULL,0,0.9)
##' normalizeVector(v,2,6,10,20)
.normalize.vector <-
  function(v,
           source.x1=NULL, source.x2=NULL,
           target.x1=0, target.x2=1, na.rm=TRUE)
{
  if (is.null(source.x1))
    source.x1 <- min(v,na.rm=na.rm)

  if (is.null(source.x2))
    source.x2 <- max(v,na.rm=na.rm)
  
  a <- (target.x2 - target.x1)/(source.x2 - source.x1)
  b <- (source.x2*target.x1 - source.x1*target.x2)/(source.x2 - source.x1)
  v <- v*a + b
  v
}


## .load.package <-
##   function(pkg,msg=NULL)
## {
##   if (is.na(packageDescription(pkg))) {
##     if (!.invalid(msg)){
##       cat(sprintf("%s\n"),msg)
##     }
##     flag <- readline(sprintf("Install \"%s\" now? [Y/n]",pkg))
##     if (tolower(flag) %in% c("","y","yes")){
##       install.packages(pkgs=pkg,...)
##     } else {
##       stop("Package \"",pkg, "\" should be installed to continue.")
##     }
##   }
##   require(pkg,quietly=TRUE,character.only=TRUE)
## }



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------



##' This function can be used to add legends to plots.  Note that a call
##' to the function \code{\link{locator}(1)} can be used in place of the \code{x}
##' and \code{y} arguments.
##' 
##' see \code{legend} in package:graphics for details;
##' Note: Old versions of graphics:::legend do not have `border' option.
##' @title Add Legends to Plots
##' @param x the x coordinates to be used to position the legend.
##' @param y the y coordinates to be used to position the legend.
##' \code{x} and \code{y} can be specified by keyword or in any way which is accepted by
##' \code{\link{xy.coords}}: See \sQuote{Details}.
##' @param legend a character or \link{expression} vector.
##' of length \eqn{\ge 1}{>= 1} to appear in the legend.  Other
##' objects will be coerced by \code{\link{as.graphicsAnnot}}.
##' @param fill if specified, this argument will cause boxes filled
##' with the specified colors (or shaded in the specified colors)
##' to appear beside the legend text.
##' @param col the color of points or lines appearing in the legend.
##' @param border the border color for the boxes (used only if \code{fill} is
##' specified).
##' @param lty the line types for lines appearing in the legend.
##' @param lwd the line widths for lines appearing in the legend.
##' One of \code{lty} and \code{lwd} \emph{must} be specified for line drawing.
##' @param pch the plotting symbols appearing in the legend, either as
##' vector of 1-character strings, or one (multi character)
##' string.  \emph{Must} be specified for symbol drawing.
##' @param angle angle of shading lines.
##' @param density the density of shading lines, if numeric and
##' positive. If \code{NULL} or negative or \code{NA} color filling
##' is assumed.
##' @param bty the type of box to be drawn around the legend.  The allowed
##' values are \code{"o"} (the default) and \code{"n"}.
##' @param bg the background color for the legend box.  (Note that this is
##' only used if \code{bty != "n"}.)
##' @param box.lwd the line type for the legend box.
##' @param box.lty the line width for the legend box.
##' @param box.col the color for the legend box.
##' @param pt.bg the background color for the \code{\link{points}},
##' corresponding to its argument \code{bg}.
##' @param cex character expansion factor \bold{relative} to current
##' \code{par("cex")}.
##' @param pt.cex expansion factor(s) for the points.
##' @param pt.lwd line width for the points, defaults to the one for
##' lines, or if that is not set, to \code{par("lwd")}.
##' @param xjust how the legend is to be justified relative to the legend
##' x location.  A value of 0 means left justified, 0.5 means centered
##' and 1 means right justified.
##' @param yjust the same as \code{xjust} for the legend y location.
##' @param x.intersp character interspacing factor for horizontal (x) spacing.
##' @param y.intersp the same for vertical (y) line distances.
##' @param adj numeric of length 1 or 2; the string adjustment for legend
##' text.  Useful for y-adjustment when \code{labels} are
##' \link{plotmath} expressions.
##' @param text.width the width of the legend text in x (\code{"user"})
##' coordinates.  (Should be positive even for a reversed x axis.)
##' Defaults to the proper value computed by \code{\link{strwidth}(legend)}.
##' @param text.col the color used for the legend text.
##' @param merge logical; if \code{TRUE}, merge points and lines but
##' not filled boxes.  Defaults to \code{TRUE} if there are points and lines.
##' @param trace logical; if \code{TRUE}, shows how \code{legend} does all
##' its magical computations.
##' @param plot logical.  If \code{FALSE}, nothing is plotted but the
##' sizes are returned.
##' @param ncol the number of columns in which to set the legend items
##' (default is 1, a vertical legend).
##' @param horiz logical; if \code{TRUE}, set the legend horizontally
##' rather than vertically (specifying \code{horiz} overrides the \code{ncol}
##' specification).
##' @param title a character string or length-one expression giving a
##' title to be placed at the top of the legend.  Other objects will be
##' coerced by \code{\link{as.graphicsAnnot}}.
##' @param inset inset distance(s) from the margins as a fraction of the
##' plot region when legend is placed by keyword.
##' @param xpd if supplied, a value of the graphical parameter 'xpd' to be
##' used while the legend is being drawn.
##' @param title.col color for \code{title}.
##' @param title.adj horizontal adjustment for \code{title}: see the help for
##' \code{par("adj")}.
##' @param seg.len the length of lines drawn to illustrate \code{lty} and/or \code{lwd}
##' (in units of character widths).
legend <- 
  function(x, y = NULL, legend, fill = NULL, col = par("col"), 
           border = "black", lty, lwd, pch, angle = 45, density = NULL, 
           bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
           box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, 
           xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0,0.5),
           text.width = NULL, text.col = par("col"), merge = do.lines && 
           has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE, 
           title = NULL, inset = 0, xpd, title.col = text.col, title.adj = 0.5, 
           seg.len = 2) 
{
  if (missing(legend) && !missing(y) && (is.character(y) || 
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  title <- as.graphicsAnnot(title)
  if (length(title) > 1) 
    stop("invalid title")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if (is.call(legend)) 
    1
  else length(legend)
  if (n.leg == 0) 
    stop("'legend' is of length 0")
  auto <- if (is.character(x)) 
    match.arg(x, c("bottomright", "bottom", "bottomleft", 
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) 
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  xlog <- par("xlog")
  ylog <- par("ylog")
  rect2 <- function(left, top, dx, dy, density = NULL, angle, 
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density, 
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    text(x, y, ...)
  }
  if (trace) 
    catn <- function(...) do.call("cat", c(lapply(list(...), 
                                                  formatC), list("\n")))
  cin <- par("cin")
  Cex <- cex * par("cex")
  if (is.null(text.width)) 
    text.width <- max(abs(strwidth(legend, units = "user", 
                                   cex = cex)))
  else if (!is.numeric(text.width) || text.width < 0) 
    stop("'text.width' must be numeric, >= 0")
  xc <- Cex * xinch(cin[1L], warn.log = FALSE)
  yc <- Cex * yinch(cin[2L], warn.log = FALSE)
  if (xc < 0) 
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
  ychar <- yextra + ymax
  if (trace) 
    catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
                                                   ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- xbox
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
                                                          0))) || !missing(lwd)
  n.legpercol <- if (horiz) {
    if (ncol != 1) 
      warning("horizontal specification overrides: Number of columns := ", 
              n.leg)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge) 
      -0.7
    else 0
  }
  else if (merge) 
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
                                                      type = "c") > 1) {
      if (length(pch) > 1) 
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
  }
  if (is.na(auto)) {
    if (xlog) 
      x <- log10(x)
    if (ylog) 
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust)) 
      xjust <- 0.5
    if (missing(yjust)) 
      yjust <- 0.5
  }
  else {
    h <- (n.legpercol + (!is.null(title))) * ychar + yc
    w0 <- text.width + (x.intersp + 1) * xchar
    if (mfill) 
      w0 <- w0 + dx.fill
    if (do.lines) 
      w0 <- w0 + (seg.len + +x.off) * xchar
    w <- ncol * w0 + 0.5 * xchar
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user", 
                                               cex = cex) + 0.5 * xchar)) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep(inset, length.out = 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = , 
                     right = usr[2L] - w - insetx, bottomleft = , 
                     left = , topleft = usr[1L] + insetx, bottom = , 
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
                    h + insety, topleft = , top = , topright = usr[4L] - 
                    insety, left = , right = , center = (usr[3L] + 
                                                         usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace) 
      catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1), 
                                              rep.int(n.legpercol, ncol)))[1L:n.leg]
  yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol, 
                                             ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
  if (mfill) {
    if (plot) {
      fill <- rep(fill, length.out = n.leg)
      rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
            col = fill, density = density, angle = angle, 
            border = border)
    }
    xt <- xt + dx.fill
  }
  if (plot && (has.pch || do.lines)) 
    col <- rep(col, length.out = n.leg)
  if (missing(lwd)) 
    lwd <- par("lwd")
  if (do.lines) {
    if (missing(lty)) 
      lty <- 1
    lty <- rep(lty, length.out = n.leg)
    lwd <- rep(lwd, length.out = n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
    if (trace) 
      catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
           yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
    if (plot) 
      segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
                col = col[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if (has.pch) {
    pch <- rep(pch, length.out = n.leg)
    pt.bg <- rep(pt.bg, length.out = n.leg)
    pt.cex <- rep(pt.cex, length.out = n.leg)
    pt.lwd <- rep(pt.lwd, length.out = n.leg)
    ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
    x1 <- (if (merge && do.lines) 
           xt - (seg.len/2) * xchar
    else xt)[ok]
    y1 <- yt[ok]
    if (trace) 
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
           ", ...)")
    if (plot) 
      points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xchar
  if (plot) {
    if (!is.null(title)) 
      text2(left + w * title.adj, top - ymax, labels = title, 
            adj = c(title.adj, 0), cex = cex, col = title.col)
    text2(xt, yt, labels = legend, adj = adj, cex = cex, 
          col = text.col)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top), 
                 text = list(x = xt, y = yt)))
}

