## lattice version of gplot.hexbin

## There are two major problems.  (1) For comparability across panels,
## we want the same mincnt and maxcnt in all panels.  However, a
## suitable default can really only be determined at printing time,
## since it would depend on the physical dimensions of the panel.  (2)
## there is no proper way to communicate the mincnt and maxcnt to the
## legend.

## Tentative solution: the counts can be calculated once enough things
## are known, namely the aspect ratio, xbins and [xy]bnds.  An
## important question then is whether [xy]bnds should be [xy]lim or
## range([xy]).  Both should be allowed, since [xy]lim makes them
## comparable, range([xy]) potentially shows more detail.  For
## relation != "same", both are more or less similar.  An important
## observation is that with range([xy]), 'shape = aspect ratio of
## panel' does not guarantee symmetric hexagons, so shape has to be
## different for each panel.

## Only feasible approach I can think of is to produce the trellis
## object first (with known aspect, so aspect="fill" is absolutely
## no-no), then analyze the limits and relevant panel arguments to get
## 'maxcnt' (essentially doing a dry run of the panel calculations).
## This needs undocumented knowledge of the trellis object, which is
## kinda not good, but at least it gets the job done.  Once we know
## maxcnt, we can also set up a suitable legend function.

## Unfortunately, this has the potential to screw up update calls that
## modify certain things.  Is there any way to capture those?  Maybe
## make a new class that inherits from "trellis".  For now, we'll
## pretend that the problem doesn't exist.


## tool borrowed from lattice
updateList <- function (x, val)
{
    if (is.null(x)) x <- list()
    modifyList(x, val)
}


prepanel.hexbinplot <-
    function(x, y, type = character(0),...)
{
   if('tmd'%in%type){
     tmp <- x
     x <- (y + x)/sqrt(2)
     y <- (y - tmp)/sqrt(2)
   }
   ans <-
     list(xlim = range(x, finite = TRUE),
          ylim = range(y, finite = TRUE),
          dx = IQR(x,na.rm=TRUE),
          dy = IQR(y,na.rm=TRUE))
}


panel.hexbinplot <-
    function(x, y, ..., groups = NULL)
{
    if (is.null(groups)) panel.hexbin(x, y, ...)
    else panel.hexpose(x, y, ..., groups = groups)
}


panel.hexbin <-
    function(x, y,
             xbins = 30,
             xbnds = c("data", "panel"), # was: xbnds = c("panel", "data"),
             ybnds = c("data", "panel"), # was: ybnds = c("panel", "data"),

             ## special args
             .prelim = FALSE,
             .cpl = current.panel.limits(),
             .xlim = .cpl$xlim,
             .ylim = .cpl$ylim,
             .aspect.ratio = 1, # default useful with splom(, panel = panel.hexbin)

             type = character(0),
             ...,
             check.erosion = FALSE)
{
    if ("tmd" %in% type) {
        tmp <- x
        x <- (y + x)/sqrt(2)
        y <- (y - tmp)/sqrt(2)
    }
    if (is.character(xbnds))
        xbnds <-
            switch(match.arg(xbnds),
                   panel = .xlim,
                   data = range(x, finite = TRUE))
    if (is.character(ybnds))
        ybnds <-
            switch(match.arg(ybnds),
                   panel = .ylim,
                   data = range(y, finite = TRUE))
    shape <-
        .aspect.ratio * (diff(ybnds) / diff(.ylim)) /
            (diff(xbnds) / diff(.xlim))
    if (!missing(check.erosion))
        warning("explicit 'check.erosion' specification ignored")
    h <- hexbin(x = x, y = y,
                xbins = xbins, shape = shape,
                xbnds = xbnds, ybnds = ybnds)
    if (.prelim)
        return(max(h@count))

    ## have to do this because grid.hexagons croaks with unrecognized
    ## arguments:
    args <- list(dat = h, check.erosion = FALSE, ...)
    keep <- names(args) %in% names(formals(grid.hexagons))

    if ('g' %in% type) panel.grid(h = -1, v = -1)
    if ('hg' %in% type) panel.hexgrid(h)

    do.call("grid.hexagons", args[keep])

    if ("r" %in% type) panel.lmline(x, y, ...)
    if ("smooth" %in% type) panel.hexloess(h,...)
    invisible()
}

panel.hexboxplot <-
    function(x, y,
             xbins = 30,
             xbnds = c("data", "panel"), # was: xbnds = c("panel", "data"),
             ybnds = c("data", "panel"), # was: ybnds = c("panel", "data"),

             ## special args
             .prelim = FALSE,
             .cpl = current.panel.limits(),
             .xlim = .cpl$xlim,
             .ylim = .cpl$ylim,
             .aspect.ratio = 1,

             type = character(0),
             cdfcut=.25,
             shadow=.05,
             ...,
             check.erosion = TRUE)
{
    if (is.character(xbnds))
        xbnds <-
            switch(match.arg(xbnds),
                   panel = .xlim,
                   data = range(x, finite = TRUE))
    if (is.character(ybnds))
        ybnds <-
            switch(match.arg(ybnds),
                   panel = .ylim,
                   data = range(y, finite = TRUE))
    shape <-
        .aspect.ratio * (diff(ybnds) / diff(.ylim)) /
            (diff(xbnds) / diff(.xlim))
    if (!missing(check.erosion))
        warning("explicit 'check.erosion' specification ignored")
    h <-hexbin(x = x, y = y,
                xbins = xbins, shape = shape,
                xbnds = xbnds, ybnds = ybnds,IDs=TRUE)

    if (.prelim)
        return(max(h@count))

    ## have to do this because grid.hexagons croaks with unrecognized
    ## arguments:
    args <- list(dat = h, check.erosion = FALSE, ...)
    keep <- names(args) %in% names(formals(grid.hexagons))
    if ('hg' %in% type) panel.hexgrid(h)
    if ('g' %in% type) panel.grid(h = -1, v = -1)
    if(shadow) {
      eh <- erode(h,cdfcut=shadow)
      h.xy <- hcell2xy(eh,check.erosion=TRUE)
      dx <- (0.5 * diff(eh@xbnds))/eh@xbins
      dy <- (0.5 * diff(eh@ybnds))/(eh@xbins * h@shape * sqrt(3))
      hexC <- hexcoords(dx, dy, sep = NULL)
      hexpolygon(h.xy$x,h.xy$y, hexC, density = density,
	       fill = NA, border = gray(.75))
    }
    eh <- erode(h,cdfcut=cdfcut)
    h.xy <- hcell2xy(eh,check.erosion=TRUE)
    dx <- (0.5 * diff(eh@xbnds))/eh@xbins
    dy <- (0.5 * diff(eh@ybnds))/(eh@xbins * h@shape * sqrt(3))
    hexC <- hexcoords(dx, dy, sep = NULL)
    hexpolygon(h.xy$x,h.xy$y, hexC, density = density,
	       fill = "green", border = gray(.75))
    med <- which.max(eh@erode)
    xnew <- h.xy$x[med]
    ynew <- h.xy$y[med]
    hexpolygon(xnew, ynew, hexC, density = density,
               fill = "red", border =gray(.25))
    invisible()
}

panel.hexpose <-
    function(x, y, groups, subscripts,
             xbins = 30,
             xbnds = c("data", "panel"), # was: xbnds = c("panel", "data"),
             ybnds = c("data", "panel"), # was: ybnds = c("panel", "data"),

             ## special args
             .prelim = FALSE,
             .cpl = current.panel.limits(),
             .xlim = .cpl$xlim,
             .ylim = .cpl$ylim,
             .aspect.ratio = 1,
             #erode Args
             cdfcut=.05,
             #hdiff Args
             hexpose.focus=c(1,2),
             hexpose.focus.colors=c("yellow","blue"),
             hexpose.focus.border=c("cyan","orange"),
             hexpose.median.color="red",
             hexpose.median.border="black",
             arrows = TRUE,
             size = unit(0.1, "inches"),
             arrow.lwd = 2,
             eps = 1e-6,
             type = character(0),
             ...,
             check.erosion = TRUE)
{
  if (is.character(xbnds))
    xbnds <-
      switch(match.arg(xbnds),
             panel = .xlim,
             data = range(x, finite = TRUE))
  if (is.character(ybnds))
    ybnds <-
      switch(match.arg(ybnds),
             panel = .ylim,
             data = range(y, finite = TRUE))
  shape <-
    .aspect.ratio * (diff(ybnds) / diff(.ylim)) /
      (diff(xbnds) / diff(.xlim))
  if (is.numeric(groups)) groups <- as.character(groups[subscripts])
  else groups <- groups[subscripts]
  binL <- hexList(x, y, given=groups, xbins=xbins, shape=shape,
                  xbnds=xbnds, ybnds=ybnds)
  if ('hs' %in% type) lapply(binL@hbins,smooth.hexbin)
  binL@hbins <- lapply(binL@hbins,erode,cdfcut=cdfcut)
  if ('hg' %in% type) panel.hexgrid(binL@hbins[[1]]) ## ???
  if ('g' %in% type) panel.grid(h = -1, v = -1)
  eroded <- unlist(lapply(binL@hbins, is, "erodebin"))
  tmph.xy <- lapply(binL@hbins, hcell2xy, check.erosion = TRUE)

  ##__________________ Construct hexagon___________________
  dx <- (0.5 * diff(binL@Xbnds))/xbins
  dy <- (0.5 * diff(binL@Ybnds))/(xbins * binL@Shape * sqrt(3))
  hexC <- hexcoords(dx = dx, dy = dy)

  ##__________________ Set up intersections and colors___________________
  ## Reorder so that the focus bin objects are at the top of the list
  if(length(hexpose.focus) < binL@n) {
    binL@hbins <- c(binL@hbins[hexpose.focus], binL@hbins[-hexpose.focus])
    binL@Bnames <- c(binL@Bnames[hexpose.focus], binL@Bnames[-hexpose.focus])
  }
  cell.stat <- all.intersect(binL@hbins)
  cell.stat.n <- apply(cell.stat, 1, sum)
  i.depth <- max(cell.stat.n)

  diff.cols <- vector(mode = "list", length = i.depth)
  levcells <- which(cell.stat.n == 1)
  whichbin <- apply(cell.stat[levcells, ], 1, which)
  ## Set all the focal colors for the unique bin cells
  ## if not specified make them equally spaced on the color wheel
  ## with high saturation and set the background bins to gray
  nfcol <- length(hexpose.focus)
  nhb <- binL@n
  nbcol <- nhb-nfcol
  fills <-
    if(is.null(hexpose.focus.colors)) {
      if(nbcol > 0)
        hsv(h = c(seq(0, 1, length = nfcol+1)[1:nfcol],rep(0, nbcol)),
            s = c(rep(1, nfcol), rep(0, nbcol)),
                        ## V = c(rep(1, nfcol), seq(.9, .1, length=nbcol))
            v = c(rep(1, nfcol), rep(.9, nbcol)))
      else hsv(h=seq(0, 1, length = nhb+1))[1:nfcol]
    }
    else {
      foc.col <- t(col2rgb(hexpose.focus.colors))/255
      if(nbcol > 0) {
        bcol <- t(col2rgb(rep(grey(.6),nbcol)))/255
        rbind(foc.col, bcol)
      }
      else foc.col
    }
  diff.cols[[1]] <- list(fill = fills, border = gray(.8))

  ##_______________ Full Cell Plotting for Unique BinL Cells_________________

    if(length(levcells) > 0) {
        for(i in unique(whichbin)) {
            pcells <-
                if(eroded[i])
                    binL@hbins[[i]]@cell[binL@hbins[[i]]@eroded]
                else binL@hbins[[i]]@cell
            pcells <- which(pcells %in% levcells[whichbin == i])

            hexpolygon(x = tmph.xy[[i]]$x[pcells],
                       y = tmph.xy[[i]]$y[pcells], hexC,
                       border = hexpose.focus.border[i] ,
                       fill = hexpose.focus.colors[i] )
        }
    }

    ## Now do the intersections. All intersections are convex
    ## combinations of the colors of the overlapping unique bins in
    ## the CIEluv colorspace.  so if the binlist is of length 2 and
    ## the focal hbins are "blue" and "yellow" respectively the
    ## intersection would be green. First I need to get this to work
    ## and then I can think about how to override this with an option
    ## in color.control. -NL

    if(i.depth > 1) {
        for(dl in 2:(i.depth)) {
            levcells <- which(cell.stat.n == dl)
            if(length(levcells) == 0) next

            whichbin <- apply(cell.stat[levcells, ], 1,
                              function(x)paste(which(x), sep = "", collapse = ":"))
            inter.nm <- unique(whichbin)
            fills <- matrix(0, length(inter.nm), 3)
            i <- 1
            for(bn in inter.nm) {
                who <- as.integer(unlist(strsplit(bn, ":")))
                ## FIXME (DS): this doesn't work
                fills[i, ] <- mixcolors2(1/length(who),
                                         diff.cols[[1]]$fill[who,])
                i <- i+1
            }
            fills <- rgb(fills[,1],fills[,2],fills[,3])
            diff.cols[[dl]] <- list(fill = fills,
                                    border = gray((i.depth-dl)/i.depth))
            ##____Full Cell Plotting for Intersecting Cells at Intersection Depth i____
            i <- 1
            for(ints in inter.nm) {
                bin.i <- as.integer(unlist(strsplit(ints, ":"))[1])
                pcells <-
                    if(eroded[bin.i])
                        binL@hbins[[bin.i]]@cell[binL@hbins[[bin.i]]@eroded]
                    else binL@hbins[[bin.i]]@cell
                pcells <- which(pcells %in% levcells[whichbin == ints])
                hexpolygon(x = tmph.xy[[bin.i]]$x[pcells],
                           y = tmph.xy[[bin.i]]$y[pcells], hexC,
                           border = diff.cols[[dl]]$border ,
                           fill = diff.cols[[dl]]$fill[i] )
                i <- i+1
            }
        }
      }

  if(any(eroded)) {
    hmeds <- matrix(unlist(lapply(binL@hbins[eroded],
                                  function(x)unlist(getHMedian(x)))),
                    ncol = 2, byrow = TRUE)
    hexpolygon(x = hmeds[, 1], y = hmeds[, 2], hexC,
               border = hexpose.median.border,
               fill = hexpose.median.color)
    if(arrows) {
      for(i in hexpose.focus) {
        for(j in hexpose.focus[hexpose.focus < i]) {
          if(abs(hmeds[i, 1] - hmeds[j, 1]) +
             abs(hmeds[i, 2] - hmeds[j, 2]) > eps)
            grid.arrows(c(hmeds[i, 1], hmeds[j, 1]),
                        c(hmeds[i, 2], hmeds[j, 2]),
                        default.units = "native",
                        length = size, gp = gpar(lwd = arrow.lwd))
        }
      }
    }
  }
  invisible()
}


hexbinplot <- function(x, data, ...) UseMethod("hexbinplot")


hexbinplot.formula <-
    function(x, data = NULL,
             prepanel = prepanel.hexbinplot,
             panel = panel.hexbinplot,
             groups = NULL,
             aspect = "xy",
             trans = NULL,
             inv = NULL,
             colorkey = TRUE,
             ...,
             maxcnt,
             legend = NULL,
             legend.width = TRUE,
             subset = TRUE)
{
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(hexbinplot)
    ccall <- match.call()
    if (is.logical(legend.width)) legend.width <- 1.2 * as.numeric(legend.width)
    if (is.character(aspect) && aspect == "fill")
        stop("aspect = 'fill' not permitted")
    if (!is.null(trans) && is.null(inv))
        stop("Must supply the inverse transformation 'inv'")
    ccall$data <- data
    ccall$prepanel <- prepanel
    ccall$panel <- panel
    ccall$aspect <- aspect
    ccall$trans <- trans
    ccall$inv <- inv
    ccall$legend <- legend
    ccall[[1]] <- quote(lattice::xyplot)
    ans <- eval(ccall, parent.frame())

    ## panel needs to know aspect ratio to calculate shape
    ans <- update(ans, .aspect.ratio = ans$aspect.ratio)

    ## also need maxcnt, o.w. can't draw legend, panels not comparable
    ## either
    if (missing(maxcnt))
        maxcnt <-
            max(mapply(panel.hexbinplot, ## note: not 'panel'
                       x = lapply(ans$panel.args, "[[", "x"),
                       y = lapply(ans$panel.args, "[[", "y"),
                       .xlim =
                       if (is.list(ans$x.limits)) ans$x.limits
                       else rep(list(ans$x.limits), length(ans$panel.args)),
                       .ylim =
                       if (is.list(ans$y.limits)) ans$y.limits
                       else rep(list(ans$y.limits), length(ans$panel.args)),
                       MoreArgs =
                       c(ans$panel.args.common,
                         list(.prelim = TRUE, .cpl = NA))))
    ans <- update(ans, maxcnt = maxcnt)
    if (colorkey)
        ans <-
            update(ans,
                   legend = updateList(ans$legend,
                   list(right =
                        list(fun = hexlegendGrob,
                             args =
                             list(maxcnt = maxcnt,
                                  trans = trans,
                                  inv = inv,
                                  legend = legend.width,
                                  ...)))))
    ans$call <- ocall
    ans
}



old.hexbinplot.formula <-
    function(x, data = parent.frame(),
             prepanel = prepanel.hexbinplot,
             panel = if (is.null(groups)) panel.hexbinplot
                    else panel.hexpose,
             groups=NULL,
             aspect = "xy",
             trans = NULL,
             inv = NULL,
             colorkey = TRUE,
             ...,
             maxcnt,
             legend = NULL,
             legend.width = TRUE)
{
    if (is.logical(legend.width))
        legend.width <- 1.2 * as.numeric(legend.width)
    if (is.character(aspect) && aspect == "fill")
        stop("aspect = 'fill' not permitted")
    if (!is.null(trans) && is.null(inv))
        stop("Must supply the inverse transformation 'inv'")
    groups <- eval(substitute(groups), data, parent.frame())
    ## There must be a better way to handle this, ugh.
    ans <-
        if(is.null(groups))
        {
            xyplot(x, data = data,
                   prepanel = prepanel,
                   panel = panel,
                   aspect = aspect,
                   trans = trans,
                   inv = inv,
                   legend = legend,
                   ...)
        }
        else
        {
            xyplot(x, data = data,
                   prepanel = prepanel,
                   panel = panel,
                   groups=groups,
                   aspect = aspect,
                   trans = trans,
                   inv = inv,
                   legend = legend,
                   ...)
        }
    ## panel needs to know aspect ratio to calculate shape
    ans <- update(ans, .aspect.ratio = ans$aspect.ratio)

    ## also need maxcnt, o.w. can't draw legend, panels not comparable
    ## either
    if (missing(maxcnt))
        maxcnt <-
            max(mapply(panel.hexbinplot, ## note: not 'panel'
                       x = lapply(ans$panel.args, "[[", "x"),
                       y = lapply(ans$panel.args, "[[", "y"),
                       .xlim =
                       if (is.list(ans$x.limits)) ans$x.limits
                       else rep(list(ans$x.limits), length(ans$panel.args)),
                       .ylim =
                       if (is.list(ans$y.limits)) ans$y.limits
                       else rep(list(ans$y.limits), length(ans$panel.args)),
                       MoreArgs =
                       c(ans$panel.args.common,
                         list(.prelim = TRUE, .cpl = NA))))
    ans <- update(ans, maxcnt = maxcnt)
    if (colorkey)
        ans <-
            update(ans,
                   legend = updateList(ans$legend,
                   list(right =
                        list(fun = hexlegendGrob,
                             args =
                             list(maxcnt = maxcnt,
                                  trans = trans,
                                  inv = inv,
                                  legend = legend.width,
                                  ...)))))
    ans
}


## want a grob instead of actual plotting

hexlegendGrob <-
    function(legend = 1.2,
             inner = legend / 5,
             cex.labels = 1,
             cex.title = 1.2,
             style = "colorscale",
             minarea = 0.05, maxarea = 0.8,
             mincnt = 1, maxcnt,
             trans = NULL, inv = NULL,
             colorcut = seq(0, 1, length = 17),
             density = NULL, border = NULL, pen = NULL,
             colramp = function(n) { LinGray(n,beg = 90,end = 15) },
             ...,
             vp = NULL,
             draw = FALSE)
{
    ## the formal arg matching should happen
    style <- match.arg(style, eval(formals(grid.hexagons)[["style"]]))
    if (style %in% c("centroids", "lattice", "colorscale")) {
        ## _______________tranformations_______________________
        if(is.null(trans))
        {
            sc <- maxcnt - mincnt
            bnds <- round(mincnt + sc * colorcut)
        }
        else
        {
            if(!is.function(trans) && !is.function(inv))
                stop("'trans' and 'inv' must both be functions if 'trans' is not NULL")
            con <- trans(mincnt)
            sc <- trans(maxcnt) - con
            bnds <- round(inv(con + sc * colorcut))
        }
    }

    ## grob
    ans <-
        switch(style,
               "colorscale" = {

                   n <- length(bnds)
                   pen <- colramp(n-1)

                   ## rectangles instead of polygons
                   ## pol <-
                   ##     rectGrob(x = 0.5, y = 1:(n-1)/n,
                   ##              height = 1/n,
                   ##              default.units = "npc",
                   ##              gp = gpar(fill = pen, col = border))

                   hexxy <- hexcoords(dx = 1, n = 1)[c("x", "y")]
                   maxxy <- max(abs(unlist(hexxy)))
                   hexxy <- lapply(hexxy, function(x) 0.5 * x/ maxxy)

                   pol <-
                       polygonGrob(x = 0.5 + rep(hexxy$x, n-1),
                                   y = (rep(1:(n-1), each = 6) + hexxy$y) / n,
                                   id.lengths = rep(6, n-1),
                                   gp = gpar(fill = pen, col = border),
                                   default.units = "npc")
                   txt <-
                       textGrob(as.character(bnds),
                                x = 0.5,
                                y = (0:(n-1) + 0.5) / n,
                                gp = gpar(cex = cex.labels),
                                default.units = "npc")
                   ttl <- textGrob("Counts", gp = gpar(cex = cex.title))

                   key.layout <-
                       grid.layout(nrow = 2, ncol = 2,
                                   heights =
                                   unit(c(1.5, 1),
                                        c("grobheight", "grobheight"),
                                        data = list(ttl, txt)),
                                   widths =
                                   unit(c(1/n, 1),
                                        c("grobheight", "grobwidth"),
                                        data = list(pol, txt)),
                                   respect = TRUE)
                   key.gf <- frameGrob(layout = key.layout, vp = vp)
                   key.gf <- placeGrob(key.gf, ttl, row = 1, col = 1:2)
                   key.gf <- placeGrob(key.gf, pol, row = 2, col = 1)
                   key.gf <- placeGrob(key.gf, txt, row = 2, col = 2)
                   key.gf
               },
               "centroids" = ,
               "lattice" = {
                   warning("legend shows relative sizes")

                   ## Note: it may not be impossible to get absolute
                   ## sizes.  The bigger problem is that when
                   ## [xy]bnds="data", the sizes (for the same count) may
                   ## not be the same across panels.  IMO, that's a more
                   ## useful feature than getting the absolute sizes
                   ## right.

                   radius <- sqrt(minarea + (maxarea - minarea) * colorcut)
                   n <- length(radius)
                   if(is.null(pen)) pen <- 1
                   if(is.null(border)) border <- pen

                   hexxy <- hexcoords(dx = 1, n = 1)[c("x", "y")]
                   maxxy <- max(abs(unlist(hexxy)))
                   hexxy <- lapply(hexxy, function(x) 0.5 * x/ maxxy)

                   pol <-
                       polygonGrob(x = 0.5 + rep(radius, each = 6) * rep(hexxy$x, n),
                                   y = (rep(0.5 + 1:n, each = 6) +
                                        rep(radius, each = 6) * hexxy$y - 1) / n,
                                   id.lengths = rep(6, n),
                                   gp = gpar(fill = pen, col = border),
                                   default.units = "npc")
                   txt <-
                       textGrob(as.character(bnds),
                                x = 0.5,
                                y = (1:n - 0.5) / n,
                                gp = gpar(cex = cex.labels),
                                default.units = "npc")
                   ttl <- textGrob("Counts", gp = gpar(cex = cex.title))

                   key.layout <-
                       grid.layout(nrow = 2, ncol = 2,
                                   heights =
                                   unit(c(1.5, 1),
                                        c("grobheight", "grobheight"),
                                        data = list(ttl, txt)),
                                   widths =
                                   unit(c(1/n, 1),
                                        c("grobheight", "grobwidth"),
                                        data = list(pol, txt)),
                                   respect = TRUE)
                   key.gf <- frameGrob(layout = key.layout, vp = vp)

                   key.gf <- placeGrob(key.gf, ttl, row = 1, col = 1:2)
                   key.gf <- placeGrob(key.gf, pol, row = 2, col = 1)
                   key.gf <- placeGrob(key.gf, txt, row = 2, col = 2)
                   key.gf
               },
               "nested.lattice" = ,
               "nested.centroids" = {
                   dx <- inner/2
                   dy <- dx/sqrt(3)
                   hexC <- hexcoords(dx, dy, n = 1, sep = NULL)

                   ## _____________x scaling_____________________________
                   numb <- cut(floor(legend/inner), breaks = c(-1, 0, 2,4))
                   ## Note: In old code
                   ##	top breaks=c(-1,0,2,4,8), numb<- 5 and size=1:9
                   if (is.na(numb)) numb <- 4
                   switch(numb,
                      {
                          warning("not enough space for legend")
                          return(textGrob(""))
                      },
                          size <- 5,
                          size <- c(1, 5, 9),
                          size <- c(1, 3, 5, 7, 9))
                   xmax <- length(size)
                   radius <- sqrt(minarea + (maxarea - minarea) * (size - 1)/9)
                   txt <- as.character(size)
                   ##___________________y scaling_____________________
                   lab <- c("Ones", "Tens", "Hundreds",
                            "Thousands", "10 Thousands", "100 Thousands",
                            "Millions", "10 Millions",
                            "100 Millions", "Billions")
                   power <- floor(log10(maxcnt)) + 1
                   yinc <- 16 * dy
                   ysize <- yinc * power
                   xmid <- 0
                   x <- inner * (1:xmax - (1 + xmax)/2) + xmid
                   n <- length(x)
                   tx <- rep.int(hexC$x, n)
                   ty <- rep.int(hexC$y, n)
                   six <- rep.int(6:6, n)
                   ## y <- rep.int(3 * dy - yinc, xmax)
                   y <- rep.int(3 * dy - 0.75 * yinc, xmax)

                   if (is.null(pen)) {
                       pen <- 1:power +1
                       pen <- cbind(pen, pen +10)
                   }
                   if (is.null(border)) border <- TRUE

                   key.layout <-
                       grid.layout(nrow = 1, ncol = 1,
                                   heights = unit(ysize, "inches"),
                                   widths = unit(legend, "inches"),
                                   respect = TRUE)
                   key.gf <- frameGrob(layout = key.layout, vp = vp)

                   ## for debugging
                   ## key.gf <-
                   ##     placeGrob(key.gf, rectGrob(gp = gpar(fill = "transparent")))

                   n6 <- rep.int(6, n)
                   for(i in 1:power) {
                       y <- y + yinc
                       key.gf <-
                           placeGrob(key.gf,
                                     polygonGrob(x = unit(legend / 2 + rep.int(hexC$x, n) + rep.int(x, n6), "inches"),
                                                 y = unit(rep.int(hexC$y, n) + rep.int(y, n6), "inches"),
                                                 id.lengths = n6,
                                                 gp =
                                                 gpar(col = pen[i, 1],
                                                      fill = if (border) 1 else pen[i, 1])),
                                     row = 1, col = 1)

                       key.gf <-
                           placeGrob(key.gf,
                                     polygonGrob(x = legend / 2 + tx * rep.int(radius, six) + rep.int(x, six),
                                                 y = ty * rep.int(radius, six) + rep.int(y, six),
                                                 default.units = "inches", id=NULL,
                                                 id.lengths=rep(6,n),
                                                 gp = gpar(fill = pen[i,2], col = border)),
                                     row = 1, col = 1)

                       key.gf <-
                           placeGrob(key.gf,
                                     textGrob(txt,
                                              x = legend / 2 + x,
                                              y = y - 4.5 * dy,
                                              default.units = "inches",
                                              gp = gpar(cex = cex.labels)),
                                     row = 1, col = 1)
                       key.gf <-
                           placeGrob(key.gf,
                                     textGrob(lab[i],
                                              x = legend / 2 + xmid,
                                              y = y[1] + 4.5 * dy,
                                              default.units = "inches",
                                              gp = gpar(cex = 1.3 * cex.title)),
                                     row = 1, col = 1)
                   }
                   key.gf
               })
    if (draw)
    {
        grid.draw(ans)
        invisible(ans)
    }
    else ans
}
