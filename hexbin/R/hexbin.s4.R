## namespace *internal* function:
addBit <- function(bnds, f = 0.05) bnds + c(-f, f) * diff(bnds)
hexbin <-
    function(x, y = NULL, xbins = 30, shape = 1,
	     xbnds = range(x), ybnds = range(y),
	     xlab = NULL, ylab = NULL, IDs = FALSE)
{
    call <- match.call()
    ## (x,y, xlab, ylab) dealing
    xl <- if (!missing(x)) deparse(substitute(x))
    yl <- if (!missing(y)) deparse(substitute(y))
    xy <- xy.coords(x, y, xl, yl)
    ch0 <- function(u) if(is.null(u)) "" else u
    xlab <- if (is.null(xlab)) ch0(xy$xlab) else xlab
    ylab <- if (is.null(ylab)) ch0(xy$ylab) else ylab
    if(! (is.character(xlab) || is.expression(xlab)))
        stop("xlab must be a character or expression")
    if(! (is.character(ylab) || is.expression(ylab)))
        stop("ylab must be a character or expression")

    x <- xy$x
    y <- xy$y
    n <- length(x)
    na <- is.na(x) | is.na(y)
    has.na <- any(na)
    if (has.na) {
	ok <- !na
	x <- x[ok]
	y <- y[ok]
        n0 <- n
        na.pos <- which(na)
	n <- length(x)
    }
    if(diff(xbnds) <= 0)
	stop("xbnds[1] < xbnds[2] is not fulfilled")
    if(!missing(xbnds) && any(sign(xbnds - range(x)) == c(1,-1)))
	stop("'xbnds' must encompass range(x)")
    if(diff(ybnds) <= 0)
	stop("ybnds[1] < ybnds[2] is not fulfilled")
    if(!missing(ybnds) && any(sign(ybnds - range(y)) == c(1,-1)))
	stop("'ybnds' must encompass range(y)")
    jmax <- floor(xbins + 1.5001)
    #imax <- 2 * floor((xbins * shape)/sqrt(3) + 1.5001)
    c1 <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
    imax <- trunc((jmax*c1 -1)/jmax + 1)
    lmax <- jmax * imax
    ans <- .Fortran("hbin",
	      x = as.double(x),
	      y = as.double(y),
	      cell = integer(lmax),
	      cnt = integer(lmax),
	      xcm = double(lmax),
	      ycm = double(lmax),
	      xbins = as.double(xbins),
	      shape = as.double(shape),
	      xbnds = as.double(xbnds),
	      ybnds = as.double(ybnds),
	      dim = as.integer(c(imax, jmax)),
	      n = as.integer(n),
	      cID = if(IDs) integer(n) else as.integer(-1),
	      PACKAGE = "hexbin")[-(1:2)]

    ## cut off extraneous stuff
    if(!IDs) ans$cID <- NULL
    if(IDs && has.na) {
      ok <- as.integer(ok)
      ok[!na] <- ans$cID
      ok[na] <- NA
      ans$cID <- ok
    }
    nc <- ans$n
    length(ans$cell) <- nc
    length(ans$cnt) <- nc
    length(ans$xcm) <- nc
    length(ans$ycm) <- nc
    if(sum(ans$cnt) != n) warning("Lost counts in binning")
    new("hexbin",
	cell = ans$cell, count = ans$cnt,
	xcm = ans$xcm, ycm = ans$ycm, xbins = ans$xbins,
	shape = ans$shape, xbnds = ans$xbnds , ybnds = ans$ybnds,
	dimen = c(imax, jmax), n = n, ncells = ans$n,
	call = call, xlab = xlab, ylab = ylab, cID = ans$cID, cAtt = integer(0))
    #dimen = ans$dim
}## hexbin

setClassUnion("integer or NULL",# < virtual class, used in 'cID' slot
	      members = c("integer","NULL"))
## MM: I've learned that we should think twice before defining such
##     "or NULL" classes:
## setClassUnion("vector or NULL",# < virtual class, used in 'cAtt' slot
## 	      members = c("vector","NULL"))

setClass("hexbin",
	 representation(cell = "integer", count = "numeric",##count = "integer",
			xcm = "numeric", ycm = "numeric", xbins = "numeric",
			shape = "numeric", xbnds = "numeric",
			ybnds = "numeric", dimen = "numeric",
			n = "integer", ncells = "integer", call = "call",
                        xlab = "vector", ylab = "vector",
			#xlab = "character", ylab = "character",
			cID = "integer or NULL", cAtt = "vector")## "or NULL"
	 )


#setIs("hexbin", function(hbin) class(hbin)=="hexbin")

## FIXME: add 'validity checking method!

setGeneric("hcell2xy", function(hbin, check.erosion = TRUE)
           standardGeneric("hcell2xy"))
setMethod("hcell2xy", "hexbin", function(hbin, check.erosion = TRUE)
{
    xbins <- hbin@xbins
    xbnds <- hbin@xbnds
    c3 <- diff(xbnds)/xbins
    ybnds <- hbin@ybnds
    c4 <- (diff(ybnds) * sqrt(3))/(2 * hbin@shape * xbins)
    jmax <- hbin@dimen[2]
    cell <- hbin@cell - 1
    i <- cell %/% jmax
    j <- cell %% jmax
    y <- c4 * i + ybnds[1]
    x <- c3 * ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
    if(check.erosion && inherits(hbin,"erodebin"))
      list(x = x[hbin@eroded], y = y[hbin@eroded])
    else
      list(x = x, y = y)
})

setGeneric("getHexDxy", function(hbin) standardGeneric("getHexDxy"))
setMethod("getHexDxy", "hexbin", function(hbin){
    sx <- hbin@xbins/diff(hbin@xbnds)
    sy <- (hbin@xbins * hbin@shape)/diff(hbin@ybnds)
    list(dx=.5/sx, dy=(1/sqrt(3))/(2*sy))
})


setClass("erodebin", representation("hexbin",
                                    eroded = "logical",
                                    cdfcut = "numeric",
                                    erode = "integer"))

setGeneric("erode", function(hbin, cdfcut = 0.5) standardGeneric("erode"))

## currently define the 'hexbin' method (also) as standalone function:
erode.hexbin <- function(hbin, cdfcut = 0.5)
{
    if(!is(hbin,"hexbin")) stop("first argument must be a hexbin object")
    #bin.att <- attributes(hbin)
    cell <- hbin@cell
    cnt <- hbin@count
    tmp <- sort(cnt)
    cdf <- cumsum(tmp)/sum(cnt)
    good <- cdfcut <= cdf
    if(!any(good))
	return("no cells selected")
    crit <- min(tmp[good])
    good <- crit <= cnt
    cell <- cell[good]
    cnt <- cnt[good]
    #hbin@cell <- cell
    #hbin@count <- cnt
    n <- length(cell)
    bdim <- hbin@dimen
    L <- bdim[1] * bdim[2]
    ans <- .Fortran("herode",
		    cell  = as.integer(cell),
		    cnt	  = as.integer(cnt),
		    n	  = n,
		    bdim  = as.integer(bdim),
		    erode = integer(L),
		    ncnt  = integer(L),
		    ncell = integer(L),
		    sides = integer(L),
		    neib  = integer(6 * L),
		    exist = logical(L + 1),
		    PACKAGE = "hexbin") $ erode
    length(ans) <- n
    ehbin <- new("erodebin", hbin, cdfcut = cdfcut, eroded = good, erode = ans)
    #hbin@erode <- ans
    #class(hbin) <- c(class(hbin),"erodebin")
    ehbin
}
setMethod("erode", "hexbin", erode.hexbin)

setGeneric("getHMedian", function(ebin) standardGeneric("getHMedian"))
setMethod("getHMedian", "erodebin", function(ebin)
      {
          xy <- hcell2xy(ebin)
          stopifnot(1 == length(med <- which.max(ebin@erode)))
          med.x <- xy$x[med]
          med.y <- xy$y[med]

          list(x = med.x, y = med.y)
      })

## Still define the 'hexbin' plot method (also) as standalone function:
## This is deprecated!
gplot.hexbin <-
    function(x, style = "colorscale",
	     legend = 1.2, lcex = 1,
	     minarea = 0.04, maxarea = 0.8, mincnt = 1, maxcnt = max(x@count),
	     trans = NULL, inv = NULL,
	     colorcut = seq(0, 1, length = min(17, maxcnt)),
	     border = NULL, density = NULL, pen = NULL,
	     colramp = function(n) LinGray(n, beg = 90, end = 15),
	     xlab = NULL, ylab = NULL, main = "", newpage = TRUE,
	     type = c("p", "l", "n"), xaxt = c("s", "n"), yaxt = c("s", "n"),
	     clip="on", verbose = getOption("verbose"))
{
    if(!is(x,"hexbin"))
	stop("first argument must be a hexbin object")
    if(minarea < 0)
	stop("Minimum area must be non-negative")
    if(maxarea > 1)
	warning("Maximum area should be <= 1 this leads to overlapping hexagons")
    if(minarea > maxarea)
	stop("Minarea must be <= maxarea")
    if (length(colorcut) > 1) { # a sequence 0,...,1
	if(colorcut[1] != 0)
	    stop("Colorcut lower boundary must be 0")
	if(colorcut[length(colorcut)] != 1)
	    stop("Colorcut upper boundary must be 1")
    }
    else {
	colorcut <-
	    if(colorcut > 1) seq(0, 1, length = min(c(17, colorcut, maxcnt)))
	    else 1
    }

    if(is.logical(legend)) {
	if(legend)
	    stop("Give the legend width")
	else legend <- 0
    } else stopifnot(is.numeric(legend) && length(legend) == 1)

    type <- match.arg(type)
    xaxt <- match.arg(xaxt)
    yaxt <- match.arg(yaxt)

    ## ----- plotting starts ------------------------
    if (newpage) grid.newpage()
    hv.ob <- hexViewport(x, offset = unit(legend,"inches"))
    pushViewport(hv.ob@hexVp.off)
    grid.rect()
    if(xaxt != "n") grid.xaxis()
    if(yaxt != "n") grid.yaxis()
    ## xlab, ylab, main :
    if(is.null(xlab)) xlab <- x@xlab
    if(is.null(ylab)) ylab <- x@ylab
    if(nchar(xlab) > 0)
      grid.text(xlab, y = unit(-2, "lines"), gp = gpar(fontsize = 16))
    if(nchar(ylab) > 0)
      grid.text(ylab, x = unit(-2, "lines"), gp = gpar(fontsize = 16), rot = 90)
    if(nchar(main) > 0)
      grid.text(main, y = unit(1, "npc") + unit(1.5, "lines"),
		gp = gpar(fontsize = 18))
    if(type != "n") {
        if(clip == "on") {
            popViewport()
            pushViewport(hv.ob@hexVp.on)
        }
        grid.hexagons(x, style = style, minarea = minarea, maxarea = maxarea,
		      mincnt = mincnt, maxcnt = maxcnt, check.erosion = FALSE,
		      trans = trans, colorcut = colorcut, density = density,
		      border = border, pen = pen,
		      colramp = colramp, verbose = verbose)
    }

    popViewport()# plot
    #popViewport()# fig
    ## ----- Legend ------------------------
    if(legend > 0) {
	if(!is.null(trans) && is.null(inv))
	    stop("Must supply the inverse transformation")
	if(verbose)
	    cat("plot.hexbin( legend > 0):  ... hex.legend()\n")
        inner <- getPlt(hv.ob, ret.unit = "inches", numeric = TRUE)[1]/x@xbins
        ##inner <- as.numeric(convertUnit(hv.ob@plt[1],"inches"))/x@xbins
	##outer <- (inner * sqrt(3))/2
	##switch(style,
	##	 lattice = ,
	##	 centroids = {
	##	     if(length(colorcut) * outer > ysize - 1) {
	##		 warning("Colorcut is being shortened")
	##		 colorcut <- seq(0, 1,
	##				 max(1, floor((ysize - 1)/outer)))
	##	     }
	##	 }
	##	 )
        ysize <- getPlt(hv.ob, ret.unit = "inches", numeric = TRUE)[2]
        #as.numeric(convertUnit(hv.ob@plt[2],"inches"))
	legVp <- viewport(x = unit(1,"npc") -
                          convertX(unit(legend,"inches"), "npc"),
			  #y = convertY(unit(mai[1],"inches"),"npc"),
                          y = hv.ob@mar[1],
			  #height = unit(1,"npc") -
			      #convertY(unit(mai[3]+mai[1],"inches"),"npc"),
                          height = unit(1,"npc")-(hv.ob@mar[1]+ hv.ob@mar[3]),
			  width = convertUnit(unit(legend,"inches"),"npc"),
			  default.units = "native",
			  just = c("left","bottom"),
			  xscale = c(0, legend),
			  yscale = c(0, ysize))
	if(type != "n") {
	    pushViewport(legVp)
	    grid.hexlegend(legend, ysize = ysize, lcex = lcex, inner = inner,
			   style = style, minarea = minarea, maxarea = maxarea,
			   mincnt = mincnt, maxcnt = maxcnt,
                           trans = trans, inv = inv, colorcut = colorcut,
			   density = density, border = border, pen = pen,
			   colramp = colramp)
	    popViewport()
	}
    }

    invisible(list(plot.vp = hv.ob, legend.vp = if(legend) legVp))
} ## gplot.hexbin()

setMethod("plot", signature(x = "hexbin", y = "missing"), gplot.hexbin)

setMethod("show", "hexbin",
	  function(object) {
	      cat("'hexbin' object from call:", deparse(object@call), "\n")
	      dm <- object@dimen
	      cat("n =", object@n, " points in	nc =", object@ncells,
		  " hexagon cells in grid dimensions ", dm[1], "by", dm[2],"\n")
	      invisible(object)
	  })

setMethod("summary", "hexbin",
	  function(object, ...) {
	      show(object, ...)
	      print(summary(data.frame(cell = object@cell, count = object@count,
				       xcm = object@xcm, ycm = object@ycm),
			    ...))
	      if(!is.null(object@cID)) {
		  cat("IDs: "); str(object@cID)
	      }
	  })



if(FALSE) { ##-- todo --
#setMethod("identify"
identify.hexbin <- function(x, labels = x$cnt, offset = 0, ...)
{
    if(length(labels) != x$n)
	stop("labels not the same length as number of cells")
    ##NL: Should this be a warning?

    ## -> typically default method:
    identify(hcell2xy(x), labels = labels, offset = offset, ...)
}
}#not yet
