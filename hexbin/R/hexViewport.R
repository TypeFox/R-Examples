setOldClass("unit")
setOldClass("viewport")

smartBnds <- function(hbin, eps=.05)
{
  hxy <- hcell2xy(hbin)
  xr <- range(hxy$x)
  yr <- range(hxy$y)
  dx <- diff(xr)
  dy <- diff(yr)
  lambda <- function(a) pmax(log(a), 1)
  epsx <- c(-1,1)*(dx*eps/lambda(dx))
  epsy <- c(-1,1)*(dy*eps/lambda(dy))
  sx <- hbin@xbins/diff(hbin@xbnds)
  sy <- (hbin@xbins * hbin@shape)/diff(hbin@ybnds)
  inner <- 0.5
  outer <- 1/sqrt(3)
  dx <- inner/sx
  dy <- outer/sy
  #xb <- dx/(hbin@xbins+1)
  #yb <- dy/((1/sqrt(3))*(hbin@xbins+1)*hbin@shape)
  list(xr = xr+ c(-dx,dx)+ epsx,
       yr = yr+ c(-dy,dy)+ epsy)
}

rname <- function(n, chars = letters)
{
	## random name with  n  characters
    paste(sample(chars, size = n, replace = TRUE), collapse="")
}

setClass("hexVP",
         representation(hexVp.on = "viewport", hexVp.off = "viewport",
                        mar = "unit", fig = "unit", plt = "unit",
                        xscale = "numeric", yscale = "numeric",shape="numeric",
                        hp.name="character")
         )

hexViewport <-
function(x, offset = unit(0,"inches"), mar = NULL,
	     xbnds = NULL, ybnds = NULL, newpage = FALSE,
         clip ="off", vp.name=NULL)
{
    if(!is(x,"hexbin"))
		stop("first argument must be a hexbin object.")
    stopifnot(is.unit(offset))

    hvp <- new("hexVP")
    if (newpage)
		grid.newpage()

    if(is.null(mar)) {
		mar <- unit(0.1 + c(5,4,4,2),"lines")
    }
    else {
		if(!is.unit(mar)) stop("'mar' must be specified in unit()s")
		if(length(mar) == 1)
			mar <- rep(mar, 4)
		else if(length(mar) != 4)
			stop("'mar' must have length 1 or 4")
    }
    ## in both cases
    mai <- as.numeric(convertUnit(mar, "inches"))
    vpin <- c(convertWidth (unit(1,"npc"),"inches"), convertHeight(unit(1,"npc"),"inches"))
    fig <- c(as.numeric(convertUnit(unit(vpin[1],"inches") - offset,"inches")), as.numeric(vpin[2]))
    pin <- c(fig[1]-mai[2]-mai[4], fig[2]-mai[1]-mai[3])
    xsize <- pin[1]
    ysize <- pin[2]

    ## The point is to optimize the placement
    ## and plotting area of the plotting window with
    ## the constraint that the margins are preserved
    ## to within some epsilon. This is going to get even
    ## harder for cases where the complex layouts are
    ## being constructed. NL -- I think it is fixed now (NL --3/22/2005)

    ## Now find the maximum rectangle in fig that
    ## has the correct aspect ratio and does not spill over epsilon into
    ## the margins, i.e.  ysize/xsize - aspect.ratio < eps and
    ##                    xsize < fig[1],  ysize < fig[2]

    if(x@shape * xsize <= ysize) {
		##center <- (ysize - x@shape * xsize)/2
		center <- (ysize - x@shape * xsize)/2
		mai[1] <- mai[1] + center
		mai[3] <- mai[3] + center
		ysize <- x@shape * xsize
    } else {
		center <- (xsize - ysize/x@shape)/2
		mai[2] <- mai[2] + center
		mai[4] <- mai[4] + center
		xsize <- ysize/x@shape
    }
    ##fig <- c(pin[1]+mai[2]+ mai[4],fig[2])
    pin <- c(xsize,ysize)
    mar <- c(convertUnit(unit(mai[1],"inches"),"lines"),
	     convertUnit(unit(mai[2],"inches"),"lines"),
	     convertUnit(unit(mai[3],"inches"),"lines"),
	     convertUnit(unit(mai[4],"inches"),"lines"))
    ##pin <- c(fig[1]-(mai[2] + mai[4]),
    ##	       fig[2]-(mai[1] + mai[3]))
    margins <- rep(as.numeric(mar), length.out = 4)
    wd <- convertUnit(unit(pin[1],"inches"),"npc")
    ## (unit(sum(margins[c(2, 4)]), "lines") +
    ##			    convertUnit(unit(legend,"inches"),"lines"))
    ## Oy, mi stupido! This is the problem, need to get the bounds right
    ## here. Fixed, do we need to guard against others stupidity and put some
    ## checks on xbnds and ybnds? (NL,4/1/2005)
    if(is.null(vp.name))
		vp.name <- rname(5)
    xyb <- smartBnds(x)
    hvp@xscale <- xs <- if(is.null(xbnds)) xyb$xr else xbnds
    hvp@yscale <- ys <- if(is.null(ybnds)) xyb$yr else ybnds
    ht <- unit(1, "npc") - unit(sum(margins[c(1,3)]), "lines")
    hvp@hexVp.off <-
        viewport(x = unit(margins[2], "lines"),
                 y = unit(margins[1], "lines"),
                 width = wd, height = ht, xscale = xs, yscale = ys,
                 just = c("left", "bottom"), default.units = "native",
                 clip = "off", name = paste(vp.name,".off",sep=""))
    hvp@hexVp.on <-
        viewport(x = unit(margins[2], "lines"),
                 y = unit(margins[1], "lines"),
                 width = wd, height = ht, xscale = xs, yscale = ys,
                 just = c("left", "bottom"), default.units = "native",
                 clip = "on", name = paste(vp.name,".on",sep=""))
    hvp@mar <- unit(mar,"lines")
    hvp@fig <- convertUnit(unit(fig,"inches"),"npc")
    hvp@plt <- convertUnit(unit(pin,"inches"),"npc")
    hvp@shape <- x@shape
    ##hvp@leg <-convertUnit(offset,"npc")
    hvp
}

## Potentially:
## setGeneric("grid:::pushViewport")
## setMethod("pushViewport", signature(x="hexVP"),
##          function(hvp) { pushViewport(hvp@hexVp) })

pushHexport <- function(hvp, clip="off")
{
    if(!is(hvp, "hexVP"))
        stop("1st argument must be 'hexVP' object")
    pushViewport(if(clip=="on") hvp@hexVp.on else hvp@hexVp.off)
}

## maybe in the future
## setMethod("push",signature("hexVP"), pushHexport)

setGeneric("getMargins", function(x, ret.unit = "npc", numeric = FALSE)
           standardGeneric("getMargins"))
setMethod("getMargins", "hexVP",
          function(x, ret.unit = "npc", numeric = FALSE){
              mar <- convertUnit(x@mar,ret.unit)
              if(numeric) as.numeric(mar) else mar
          })

setGeneric("getPlt", function(x, ret.unit = "npc", numeric = FALSE)
           standardGeneric("getPlt"))
setMethod("getPlt", "hexVP",
          function(x, ret.unit = "npc", numeric = FALSE){
              plt <- convertUnit(x@plt,ret.unit)
              if(numeric) as.numeric(plt) else plt
          })

setGeneric("getFig", function(x, ret.unit = "npc", numeric = FALSE)
           standardGeneric("getFig"))
setMethod("getFig", "hexVP",
          function(x, ret.unit = "npc", numeric = FALSE){
              fig <- convertUnit(x@fig,ret.unit)
              if(numeric) as.numeric(fig) else fig
          })

## MM doesn't think it's ok to "pollute" the generic-space
##    just for basic slot accessors :

## setGeneric("getXscale", function(x)standardGeneric("getXscale"))
## setMethod("getXscale", "hexVP", function(x){ x@xscale })

## setGeneric("getYscale", function(x)standardGeneric("getYscale"))
## setMethod("getYscale", "hexVP", function(x){ x@yscale })

hexVP.abline <- function(hvp, a = NULL, b = NULL, h = numeric(0),
                         v = numeric(0), col = 'black',
                         lty = 1, lwd = 2, ...)
{
    pushHexport(hvp, clip = 'on')
    col.line <- col
    if (!is.null(a)) {
        if (inherits(a, "lm")) {
            coeff <- coef(a)
        }
        else if (!is.null(tryCatch(coef(a), error = function(e) NULL)))
            coeff <- coef(a)
        else coeff <- c(a, b)
        if (length(coeff) == 1)
            coeff <- c(0, coeff)
        if (coeff[2] == 0)
            h <- c(h, coeff[1])
        else if (!any(is.null(coeff))) {
            xx <- current.viewport()$xscale
            yy <- current.viewport()$yscale
            x <- numeric(0)
            y <- numeric(0)
            ll <- function(i, j, k, l)
				(yy[j] - coeff[1] - coeff[2] * xx[i]) * (yy[l] - coeff[1] - coeff[2] * xx[k])
            if (ll(1, 1, 2, 1) <= 0) {
                y <- c(y, yy[1])
                x <- c(x, (yy[1] - coeff[1])/coeff[2])
            }
            if (ll(2, 1, 2, 2) <= 0) {
                x <- c(x, xx[2])
                y <- c(y, coeff[1] + coeff[2] * xx[2])
            }
            if (ll(2, 2, 1, 2) <= 0) {
                y <- c(y, yy[2])
                x <- c(x, (yy[2] - coeff[1])/coeff[2])
            }
            if (ll(1, 2, 1, 1) <= 0) {
                x <- c(x, xx[1])
                y <- c(y, coeff[1] + coeff[2] * xx[1])
            }
            if (length(x) > 0)
                grid.lines(x = x, y = y, default.units = "native",
                           gp = gpar(col = col.line, lty = lty, lwd = lwd))
        }
    }
    h <- as.numeric(h)
    v <- as.numeric(v)
    for (i in seq(along = h))
        grid.lines(y = rep(h[i], 2), default.units = "native",
                   gp = gpar(col = col.line, lty = lty, lwd = lwd))
    for (i in seq(along = v))
        grid.lines(x = rep(v[i], 2), default.units = "native",
                   gp = gpar(col = col.line, lty = lty, lwd = lwd))
    popViewport()
}

hexVP.loess <- function(hbin, hvp = NULL, span = 0.4, col = 'red', n = 200)
{
    fit <- loess(hbin@ycm ~ hbin@xcm, weights = hbin@count, span = span)
    if(!is.null(hvp)) {
        pushHexport(hvp, clip = 'on')
#        grid.lines(seq(0,16, length = n),
#                   predict(fit,seq(0,16, length = n)),
#                   gp = gpar(col = col), default.units = 'native')
 		grid.lines(seq(hbin@xbnds[1], hbin@xbnds[2], length = n),
				predict(fit,seq(hbin@xbnds[1], hbin@xbnds[2], length = n)),
				gp = gpar(col = col), default.units = 'native')
        popViewport()
    }
    invisible(fit)
}
