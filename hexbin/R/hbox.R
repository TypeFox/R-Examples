hboxplot <- function(bin, xbnds = NULL, ybnds = NULL,
		     density, border = c(0,grey(.7)),
		     pen = c(2, 3), unzoom = 1.1, clip="off", reshape = FALSE,
		     xlab = NULL, ylab = NULL, main = "")
{

    ##_______________ Collect computing constants______________

    if(!is(bin,"hexbin"))
	stop("first argument must be a hexbin object")
    h.xy <- hcell2xy(bin,check.erosion=TRUE)
    ##___zoom in scaling with expanding to avoid hexagons outside plot frame___

    if(is(bin,"erodebin")) {
	h.xy$x <- h.xy$x
	h.xy$y <- h.xy$y
	nxbnds <- if(is.null(xbnds)) range(h.xy$x) else xbnds
	nybnds <- if(is.null(ybnds)) range(h.xy$y) else ybnds
	ratiox <- diff(nxbnds)/diff(bin@xbnds)
	ratioy <- diff(nybnds)/diff(bin@ybnds)

	ratio <- max(ratioy, ratiox)
	nxbnds <- mean(nxbnds) + c(-1,1)*(unzoom * ratio * diff(bin@xbnds))/2
	nybnds <- mean(nybnds) + c(-1,1)*(unzoom * ratio * diff(bin@ybnds))/2
    }
    else {
	nxbnds <- if(is.null(xbnds)) bin@xbnds else xbnds
	nybnds <- if(is.null(ybnds)) bin@ybnds else ybnds
    }
    margins <- unit(0.1 + c(5,4,4,3),"lines")
    plot.vp <- hexViewport(bin, xbnds = nxbnds, ybnds = nybnds,
                           mar=margins, newpage = TRUE)
    pushHexport(plot.vp)
    grid.rect()
    grid.xaxis()
    grid.yaxis()
    ## xlab, ylab, main :
    if(is.null(xlab)) xlab <- bin@xlab
    if(is.null(ylab)) ylab <- bin@ylab
    if(nchar(xlab) > 0)
	grid.text(xlab, y = unit(-2, "lines"), gp= gpar(fontsize= 16))
    if(nchar(ylab) > 0)
	grid.text(ylab, x = unit(-2, "lines"), gp= gpar(fontsize= 16), rot = 90)
    if(nchar(main) > 0)
	grid.text(main, y = unit(1, "npc") + unit(1.5, "lines"),
		  gp = gpar(fontsize = 18))
    if(clip=="on") {
      popViewport()
      pushHexport(plot.vp, clip="on")
    }

    cnt <- if(is(bin,"erodebin")) bin@count[bin@eroded] else bin@count

    xbins <- bin@xbins
    shape <- bin@shape
    xnew <- h.xy$x
    ynew <- h.xy$y

    ##__________________ Construct a hexagon___________________
    dx <- (0.5 * diff(bin@xbnds))/xbins
    dy <- (0.5 * diff(bin@ybnds))/(xbins * shape * sqrt(3))
    hexC <- hexcoords(dx, dy, sep = NULL)

    ##_______________ Full Cell Plotting_____________________
    hexpolygon(xnew, ynew, hexC, density = density,
	       fill = pen[2], border = border[2])

    ##______________Plotting median___________________________

    if(!is(bin,"erodebin")) {
	## No warning here, allow non-erode above!   warning("No erode component")
    }
    else {
	med <- which.max(bin@erode)
	xnew <- xnew[med]
	ynew <- ynew[med]
	hexpolygon(xnew, ynew, hexC, density = density,
		   fill = pen[1], border = border[1])
    }
    popViewport()
    invisible(plot.vp)

}# hboxplot()
