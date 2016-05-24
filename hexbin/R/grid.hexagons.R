
hexcoords <- function(dx, dy = NULL, n = 1, sep = NULL)
{
    stopifnot(length(dx) == 1)
    if(is.null(dy)) dy <- dx/sqrt(3)
    if(is.null(sep))
        list(x = rep.int(c(dx, dx,     0, -dx, -dx,    0), n),
             y = rep.int(c(dy,-dy, -2*dy, -dy,  dy, 2*dy), n),
             no.sep = TRUE)
    else
        list(x = rep.int(c(dx, dx,     0, -dx, -dx,   0, sep), n),
             y = rep.int(c(dy,-dy, -2*dy, -dy,	dy, 2*dy, sep), n),
             no.sep = FALSE)
}

hexpolygon <-
function(x, y, hexC = hexcoords(dx, dy, n = 1), dx, dy=NULL,
	     fill = 1, border = 0, hUnit = "native", ...)
{
    ## Purpose: draw hexagon [grid.]polygon()'s  around  (x[i], y[i])_i
    ## Author: Martin Maechler, Jul 2004; Nicholas for grid

    n <- length(x)
    stopifnot(length(y) == n)
    stopifnot(is.list(hexC) && is.numeric(hexC$x) && is.numeric(hexC$y))
    if(hexC$no.sep) {
        n6 <- rep.int(6:6, n)
        if(!is.null(hUnit)) {
            grid.polygon(x = unit(rep.int(hexC$x, n) + rep.int(x, n6),hUnit),
                         y = unit(rep.int(hexC$y, n) + rep.int(y, n6),hUnit),
                         id.lengths = n6,
                         gp = gpar(col= border, fill= fill))
        }
        else {
            grid.polygon(x = rep.int(hexC$x, n) + rep.int(x, n6),
                         y = rep.int(hexC$y, n) + rep.int(y, n6),
                         id.lengths = n6,
                         gp = gpar(col= border, fill= fill))
        }
    }
    else{ ## traditional graphics polygons: must be closed explicitly (+ 1 pt)
        n7 <- rep.int(7:7, n)
        polygon(x = rep.int(hexC$x, n) + rep.int(x, n7),
                y = rep.int(hexC$y, n) + rep.int(y, n7), ...)
    }
}

grid.hexagons <-
function(dat, style = c("colorscale", "centroids", "lattice",
	     "nested.lattice", "nested.centroids", "constant.col"),
         use.count=TRUE, cell.at=NULL,
         minarea = 0.05, maxarea = 0.8, check.erosion = TRUE,
         mincnt = 1, maxcnt = max(dat@count), trans = NULL,
         colorcut = seq(0, 1, length = 17),
         density = NULL, border = NULL, pen = NULL,
         colramp = function(n){ LinGray(n,beg = 90, end = 15) },
         def.unit = "native",
         verbose = getOption("verbose"))
{
    ## Warning:	 presumes the plot has the right shape and scales
    ##		 See plot.hexbin()
    ## Arguments:
    ##	dat   = hexbin object
    ##	style = type of plotting
    ##		'centroids' =  symbol area is a function of the count,
    ##			approximate location near cell center of
    ##			mass without overplotting
    ##		'lattice' = symbol area is a function of the count,
    ##			plot at lattice points
    ##		'colorscale'   = gray scale plot,
    ##			color number determined by
    ##			transformation and colorcut,
    ##			area = full hexagons.
    ##		'nested.lattice'= plots two hexagons
    ##		   background hexagon
    ##			area=full size
    ##			color number by count in powers of 10 starting at pen 2
    ##		   foreground hexagon
    ##			area by log10(cnt)-floor(log10(cnt))
    ##			color number by count in powers of 10 starting at pen 12
    ##		'nested.centroids' = like nested.lattice
    ##			but counts < 10 are plotted
    ##
    ##	minarea =  minimum symbol area as fraction of the binning cell
    ##	maxarea =  maximum symbol area as fraction of the binning cell
    ##	mincnt	=  minimum count accepted in plot
    ##	maxcnt	=  maximum count accepted in plot
    ##	trans	=  a transformation scaling counts into [0,1] to be applied
    ##		   to the counts for options 'centroids','lattice','colorscale':
    ##			 default=(cnt-mincnt)/(maxcnt-mincnt)
    ##	colorcut=  breaks for translating values between 0 and 1 into
    ##		   color classes.  Default= seq(0,1,17),
    ##	density =  for hexagon graph paper
    ##	border	   plot the border of the hexagon, use TRUE for
    ##		   hexagon graph paper
    ## Symbol size encoding:
    ##	 Area= minarea + scaled.count*(maxarea-minarea)
    ##	 When maxarea==1 and scaled.count==1, the hexagon cell
    ##	 is completely filled.
    ##
    ##	 If small hexagons are hard to see increase minarea.
    ## For gray scale encoding
    ##	  Uses the counts scaled into [0,1]
    ##	  Default gray cutpoints seq(0,1,17) yields 16 color classes
    ##	  The color number for the first class starts at 2.
    ##	       motif coding: black 15 white puts the first of the
    ##			     color class above the background black
    ##	  The function subtracts 1.e-6 from the lower cutpoint to include
    ##	  the boundary
    ## For nested scaling see the code
    ## Count scaling alternatives
    ##
    ##	log 10 and Poisson transformations
    ##	trans <- function(cnt) log10(cnt)
    ##	   min	inv   <- function(y) 10^y
    ##
    ##	     trans <- function(cnt) sqrt(4*cnt+2)
    ##	     inv   <- function(y) (y^2-2)/4
    ## Perceptual considerations.
    ##	  Visual response to relative symbol area is not linear and varies from
    ##	  person to person.  A fractional power transformation
    ##	  to make the interpretation nearly linear for more people
    ##	  might be considered.	With areas bounded between minarea
    ##	  and 1 the situation is complicated.
    ##
    ##	  The local background influences color interpretation.
    ##	  Having defined color breaks to focus attention on
    ##	  specific countours can help.
    ##
    ## Plotting the symbols near the center of mass is not only more accurate,
    ##	  it helps to reduce the visual dominance of the lattice structure.  Of
    ##	  course higher resolution binning reduces the possible distance between
    ##	  the center of mass for a bin and the bin center.  When symbols
    ##	  nearly fill their bin, the plot appears to vibrate.  This can be
    ##	  partially controlled by reducing maxarea or by reducing
    ##	  contrast.


    ##____________________Initial checks_______________________
    if(!is(dat,"hexbin"))
		stop("first argument must be a hexbin object")
    style <- match.arg(style) # so user can abbreviate
    if(minarea <= 0)
		stop("hexagons cannot have a zero area, change minarea")
    if(maxarea > 1)
		warning("maxarea > 1, hexagons may overplot")
    ##_______________ Collect computing constants______________

    if(use.count){
      cnt <- dat@count
    }
    else{
      cnt <- cell.at
      if(is.null(cnt)){
        if(is.null(dat@cAtt)) stop("Cell attribute cAtt is null")
        else cnt <- dat@cAtt
      }
    }
    xbins <- dat@xbins
    shape <- dat@shape
    tmp <- hcell2xy(dat, check.erosion = check.erosion)
    good <- mincnt <= cnt & cnt <= maxcnt
    xnew <- tmp$x[good]
    ynew <- tmp$y[good]
    cnt <- cnt[good]
    sx <- xbins/diff(dat@xbnds)
    sy <- (xbins * shape)/diff(dat@ybnds)

	##___________Transform Counts to Radius_____________________
    switch(style,
	   "centroids" = ,
	   "lattice" = ,
	   "constant.col" =,
	   "colorscale" = {
	       if(is.null(trans)) {
                   if( min(cnt,na.rm=TRUE)< 0){
                     pcnt<- cnt + min(cnt)
                     rcnt <- {
                       if(maxcnt == mincnt) rep.int(1, length(cnt))
                       else (pcnt - mincnt)/(maxcnt - mincnt)
                     }
                   }
		   else rcnt <- {
                     if(maxcnt == mincnt) rep.int(1, length(cnt))
                     else (cnt - mincnt)/(maxcnt - mincnt)
		   }
	       }
	       else {
		   rcnt <- (trans(cnt) - trans(mincnt)) /
		       (trans(maxcnt) - trans(mincnt))
		   if(any(is.na(rcnt)))
		       stop("bad count transformation")
	       }
	       area <- minarea + rcnt * (maxarea - minarea)
	   },
	   "nested.lattice" = ,
	   "nested.centroids" = {
	       diffarea <- maxarea - minarea
	       step <- 10^floor(log10(cnt))
	       f <- (cnt/step - 1)/9
	       area <- minarea + f * diffarea
	       area <- pmax(area, minarea)
	   }
	   )
    area <- pmin(area, maxarea)
    radius <- sqrt(area)

    ##______________Set Colors_____________________________
    switch(style,
	   "centroids" = ,
       "constant.col" = ,
	   "lattice" = {
               if(length(pen)!= length(cnt)){
                   if(is.null(pen)) pen <- rep.int(1, length(cnt))
               ##else if(length(pen)== length(cnt)) break
                   else if(length(pen)== 1) pen <- rep.int(pen,length(cnt))
                   else stop("'pen' has wrong length")
               }
           },
	   "nested.lattice" = ,
	   "nested.centroids" = {
	       if(!is.null(pen) && length(dim(pen)) == 2) {
		   dp <- dim(pen)
		   lgMcnt <- ceiling(log10(max(cnt)))
		   if(dp[1] != length(cnt) && dp[1] != lgMcnt ) {
		       stop ("pen is not of right dimension")
		   }
		   if( dp[1] == lgMcnt ) {
		       ind <- ceiling(log10(dat@count)) ## DS: 'dat' was 'bin' (??)
		       ind[ind == 0] <- 1
		       pen <- pen[ind,]
		   }
		   ##else break
	       }
	       else {
		   pen <- floor(log10(cnt)) + 2
		   pen <- cbind(pen, pen+10)
	       }
	   },
	   "colorscale" = {
	       ## MM: Following is quite different from bin2d's
               nc <- length(colorcut)
	       if(colorcut[1] > colorcut[nc]){
		   colorcut[1] <- colorcut[1] + 1e-06
		   colorcut[nc] <- colorcut[nc] - 1e-06
	       } else {
		   colorcut[1] <- colorcut[1] - 1e-06
		   colorcut[nc] <- colorcut[nc] + 1e-06
	       }
	       colgrp <- cut(rcnt, colorcut,labels = FALSE)
               if(any(is.na(colgrp))) colgrp <- ifelse(is.na(colgrp),0,colgrp)
	       ##NL: colramp must be a function accepting an integer n
	       ##    and returning n colors
	       clrs <- colramp(length(colorcut) - 1)
	       pen <- clrs[colgrp]
	   }
           )

    ##__________________ Construct a hexagon___________________
    ## The inner and outer radius for hexagon in the scaled plot
    inner <- 0.5
    outer <- (2 * inner)/sqrt(3)
    ## Now construct a point up hexagon symbol in data units
    dx <- inner/sx
    dy <- outer/(2 * sy)
    rad <- sqrt(dx^2 + dy^2)
    hexC <- hexcoords(dx, dy, sep=NULL)
    ##_______________ Full Cell	 Plotting_____________________
    switch(style,
	   "constant.col" = ,
	   "colorscale" = {
	       hexpolygon(xnew, ynew, hexC,
			  density = density, fill = pen,
			  border = if(!is.null(border)) border else pen)

               ## and that's been all for these styles
	       return(invisible(paste("done", sQuote(style))))
	   },
	   "nested.lattice" = ,
	   "nested.centroids" = {
	       hexpolygon(xnew, ynew, hexC,
			  density = density,
                          fill = if (is.null(border) || border) 1 else pen[,1],
			  border = pen[,1])
	   }
	   )

    ##__________________ Symbol Center adjustments_______________
    if(style == "centroids" || style == "nested.centroids") {
	xcm <- dat@xcm[good]
	ycm <- dat@ycm[good]
	## Store 12 angles around a circle and the replicate the first
	## The actual length for these vectors is determined by using
	## factor use below
	k <- sqrt(3)/2
	cosx <- c(1, k, .5, 0, -.5, -k, -1, -k,	 -.5, 0, .5,   k, 1)/sx
	siny <- c(0, .5, k, 1,	 k, .5,	 0, -.5, -k, -1, -k, -.5, 0)/sy
	## Compute distances for differences after scaling into
	## [0,size] x [0,aspect*size]
	## Then there are size hexagons on the x axis
	dx <- sx * (xcm - xnew)
	dy <- sy * (ycm - ynew)
	dlen <- sqrt(dx^2 + dy^2)
	## Find the closest approximating direction of the 12 vectors above
	cost <- ifelse(dlen > 0, dx/dlen, 0)
	tk <- (6 * acos(cost))/pi
	tk <- round(ifelse(dy < 0, 12 - tk, tk)) + 1
	## Select the available length for the approximating vector
	hrad <- ifelse(tk %% 2 == 1, inner, outer)
	## Rad is either an inner or outer approximating radius.
	## If dlen + hrad*radius <= hrad, move the center dlen units.
	## Else move as much of dlen as possible without overplotting.
	fr <- pmin(hrad * (1 - radius), dlen) # Compute the symbol centers
	## fr is the distance for the plot [0,xbins] x [0,aspect*xbins]

	## cosx and siny give the x and y components of this distance
	## in data units
	xnew <- xnew + fr * cosx[tk]
	ynew <- ynew + fr * siny[tk]
    }
    ## ________________Sized Hexagon Plotting__________________
    ## scale the symbol by radius and add to the new center
    n <- length(radius)
    if(verbose)
	cat('length = ',length(pen),"\n", 'pen = ', pen+1,"\n")
    ##switch(style,
    ##	    centroids = ,
    ##	    lattice = {if(is.null(pen))pen <- rep.int(1, n)
    ##			else pen <- rep.int(pen, n)},
    ##	   nested.lattice = ,
    ##	    nested.centroids ={
    ##	      if(
    ##		 pen[,2] <- pen[,1] + 10
    ##	  }	      )

    ## grid.polygon() closes automatically: now '6' where we had '7':
    n6 <- rep.int(6:6, n)
    pltx <- rep.int(hexC$x, n) * rep.int(radius, n6) + rep.int(xnew, n6)
    plty <- rep.int(hexC$y, n) * rep.int(radius, n6) + rep.int(ynew, n6)
    switch(style,
	   "centroids" = ,
	   "lattice" = {
	       grid.polygon(pltx, plty, default.units=def.unit, id=NULL,
			    ## density = density,
			    id.lengths= n6,
			    gp=gpar(fill = pen, col = border))
	   },
	   "nested.lattice" = ,
	   "nested.centroids" = {
	       grid.polygon(pltx, plty, default.units=def.unit, id=NULL,
			    id.lengths= n6,
			    gp=gpar(fill = pen[,2],
			    ## density = density,
			    col=if(!is.null(border)) border else pen[,2]))

	   })

}

if(FALSE){ ## considering 'hexagons' object
    setMethod("hexagons", signature(dat="hexbin"), grid.hexagons)

    erode.hexagons <- function(ebin,pen="black",border="red"){
	print("Blank for now")
    }
}
