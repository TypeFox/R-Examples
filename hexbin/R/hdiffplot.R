
### FIXME: Need to check for bin erosion
###        or fix hcell2xy so that it checks for bin erosion.
### --- Fixed hcell2xy, probably should do the same to other accessor functions
### NL

get.xrange <- function(xy.lst, xbnds)
{
    range(unlist(lapply(xy.lst,
                        function(xy, bnd)
                        xy$x[(xy$x < max(bnd)) & (xy$x > min(bnd))],
                        xbnds)))
}

get.yrange <- function(xy.lst, ybnds)
{
    range(unlist(lapply(xy.lst,
                        function(xy, bnd)
                        xy$y[(xy$y < max(bnd)) & (xy$y > min(bnd))],
                        ybnds)))
}

make.bnds <- function(binlst, xy.lst, xbnds = NULL, ybnds = NULL)
{
    if(inherits(binlst,"hexbinList")) binlst <- binlst@hbins
    if(is.null(xbnds)) xbnds <- binlst[[1]]@xbnds
    if(is.null(ybnds)) ybnds <- binlst[[1]]@ybnds

    nxbnds <- get.xrange(xy.lst, xbnds)
    nybnds <- get.yrange(xy.lst, ybnds)

    list(xbnds = xbnds, ybnds = ybnds, nxbnds = nxbnds, nybnds = nybnds)
}

all.intersect <- function(binlist)
{
    ## This will not work if all the grids are not the same
    ## Will have to rethink this if we move to non-aligned
    ## hexagon bins. NL
    if(inherits(binlist,"hexbinList")) binlist <- binlist@hbins
    ans <- matrix(FALSE, nrow = binlist[[1]]@dimen[1]*binlist[[1]]@dimen[2],
                  ncol = length(binlist))
    for(i in 1:length(binlist)) {
        if(is(binlist[[i]], "erodebin"))
            ans[binlist[[i]]@cell[binlist[[i]]@eroded], i] <- TRUE
        else ans[binlist[[i]]@cell, i] <- TRUE
    }
    ans
}

## colordist <- function() {
## }

## MM: FIXME :  `` get(where) '' is a kludge!
# EJP: outcomment, seems obsolete?
#mixcolors <- function (alpha, color1, where = class(color1))
#{
#    alpha <- as.numeric(alpha)
#    c1 <- coords(as(color1, where))
#    na <- length(alpha)
#    n1 <- nrow(c1)
#    if(na == 1)
#        alpha <- rep(alpha, n1)
#    stopifnot(sum(alpha) == 1)
#    get(where)(t(apply(c1, 2, function(cols, alpha) alpha%*%cols, alpha)))
#
#}

mixcolors2 <- function (colors, alpha, where="hsv")
{
    # colors: an n x 3 matrix of colors
    # alpha: an n x 1 vector of color mixing coefficents
    #        sum(alpha)==1 should be a restriction?
    # where: the color space to mix in (not implemented yet)
    # The reurn value is a single hex color forming the mixture
    # This function is purely linear mixing, nolinear mixing
    # would be quite interesting since the colorspaces are not really
    # linear, ie mixing alonga manifold in LUV space.
    alpha <- as.numeric(alpha)
    na <- length(alpha)
    n1 <- nrow(colors)
    if (n1 < 2) {
     warning("need more than two colors to mix")
     colors
    }
    if(na == 1)
        alpha <- rep(alpha, n1)
    stopifnot(abs(sum(alpha)-1) <= 0.01)
    #colors <- convertColor(colors,from="sRGB",to="Lab",scale.in=1)
    mix <- t(apply(colors, 2, function(cols, alpha) alpha%*%cols, alpha))
    #convertColor(mix,from="hsv",to="hex",scale.out=1,clip=TRUE)
    hsv(mix[1],mix[2],mix[3])
}

hdiffplot <-
    function(bin1, bin2 = NULL, xbnds = NULL, ybnds = NULL,
             focus = NULL,
             col.control = list(medhex = "white", med.bord = "black",
             focus = NULL, focus.border = NULL,
             back.col = "grey"),
             arrows = TRUE, size = unit(0.1, "inches"), lwd = 2,
             eps = 1e-6, unzoom = 1.08, clip ="off", xlab = "", ylab = "",
             main = deparse(mycall), ...)
{
    ## Arguments:
    ##	bin1	: hexagon bin object or a list of bin objects
    ##  bin2	: hexagon bin object or NULL
    ##            bin objects must have the same plotting bounds and shape
    ## 	border  : plot the border of the hexagon, use TRUE for
    ##            hexagon graph paper

	## Having all the same parameters ensures that all hexbin
	## objects have the same hexagon grid, and there will be no
	## problems intersecting them. When we have a suitable solution to
	## the hexagon interpolation/intersection problem this will be relaxed.
	fixx <- xbnds
    fixy <- ybnds

    if(!inherits(bin1,"hexbinList")){
      if(is.null(bin2) & is.list(bin1)) {
        bin1 <- as(bin1,"hexbinList")
      }
      else if(is.null(bin2) & (!is.list(bin1)))
        stop(" need at least 2 hex bin objects, or a hexbinList")
      else {
        if(bin1@shape != bin2@shape)
            stop("bin objects must have same shape parameter")
        if(all(bin1@xbnds == bin2@xbnds) & all(bin1@ybnds == bin2@ybnds))
            equal.bounds <- TRUE
        else stop("Bin objects need the same xbnds and ybnds")
        if(bin1@xbins != bin2@xbins)
            stop("Bin objects need the same number of bins")
        nhb <- 2
        ## Need to make a binlist class, then can do as(bin1, bin2, "binlist")
        ## or something similar (NL)
        bin1 <- list(bin1 = bin1, bin2 = bin2)
        bin1 <- as(bin1,"hexbinList")
      }
    }
    mycall <- sys.call()
    if(length(mycall) >= 4) {
        mycall[4] <- as.call(quote(.....()))
        if(length(mycall) > 4) mycall <- mycall[1:4]
    }
    if(is.null(focus)) focus <- 1:bin1@n
    ##_______________ Collect computing constants______________
    tmph.xy <- lapply(bin1@hbins, hcell2xy, check.erosion = TRUE)

    ## Check for erode bins
    eroded <- unlist(lapply(bin1@hbins, is, "erodebin"))
    shape <- bin1@Shape
    xbins <- bin1@Xbins
    bnds <- make.bnds(bin1@hbins, tmph.xy, xbnds = fixx, ybnds = fixy)
    ratiox <- diff(bnds$nxbnds)/diff(bnds$xbnds)
    ratioy <- diff(bnds$nybnds)/diff(bnds$ybnds)
    ratio <- max(ratioy, ratiox)

    nxbnds <- mean(bnds$nxbnds) + c(-1, 1)*(unzoom * ratio * diff(bnds$xbnds))/2
    nybnds <- mean(bnds$nybnds) + c(-1, 1)*(unzoom * ratio * diff(bnds$ybnds))/2

    ##__________________ Construct plot region___________________
    hvp <- hexViewport(bin1@hbins[[1]], xbnds = nxbnds, ybnds = nybnds,
                       newpage = TRUE)
    pushHexport(hvp)
    grid.rect()
    grid.xaxis()
    grid.yaxis()
    if(nchar(xlab) > 0)
        grid.text(xlab, y = unit(-2, "lines"), gp = gpar(fontsize = 16))
    if(nchar(ylab) > 0)
        grid.text(ylab, x = unit(-2, "lines"), gp = gpar(fontsize = 16), rot = 90)
    if(sum(nchar(main)) > 0)
        grid.text(main, y = unit(1, "npc") + unit(1.5, "lines"),
                  gp = gpar(fontsize = 18))

    if(clip=='on'){
      popViewport()
      pushHexport(hvp,clip="on")
    }
    ##__________________ Construct hexagon___________________
    dx <- (0.5 * diff(bin1@Xbnds))/xbins
    dy <- (0.5 * diff(bin1@Ybnds))/(xbins * shape * sqrt(3))
    hexC <- hexcoords(dx = dx, dy = dy)

    ##__________________ Set up intersections and colors___________________
    if(length(focus) < bin1@n) {
        bin1@hbins <- c(bin1@hbins[focus], bin1@hbins[-focus])
        bin1@Bnames <- c(bin1@Bnames[focus], bin1@Bnames[-focus])
    }
    cell.stat <- all.intersect(bin1@hbins)
    cell.stat.n <- apply(cell.stat, 1, sum)
    i.depth <- max(cell.stat.n)

### I will do this as a recursive function once I get
### The colors worked out! In fact for more than three
### bin objects there is no other way to do this but recursively!!!
### NL. -- Well this solution is like recursion :)
    diff.cols <- vector(mode = "list", length = i.depth)
    levcells <- which(cell.stat.n == 1)
    whichbin <- apply(cell.stat[levcells, ], 1, which)

    ## Set all the focal colors for the unique bin cells
    ## if not specified make them equally spaced on the color wheel
    ## with high saturation and set the background bins to gray
    nfcol <- length(focus)
    nhb <- bin1@n
    nbcol <- nhb-nfcol
    fills <-
        if(is.null(col.control$focus)) {
                if(nbcol > 0)
                  matrix(c(seq(0, 360, length = nfcol+1)[1:nfcol]/360, rep(0, nbcol),
                           rep(1, nfcol), rep(0, nbcol),rep(1, nfcol), rep(.9, nbcol)),
                         ncol = 3)
                        ## V = c(rep(1, nfcol), seq(.9, .1, length=nbcol))

                else #matrix(c(seq(0, 360, length = nhb+1), s=1, v=1)[1:nfcol]
                  matrix(c(seq(0, 360, length = nhb+1)/360,
                            rep(1,nhb+1),
                            rep(1,nhb+1)), ncol = 3)[1:nhb,]
        }
        else {
            foc.col <- t(rgb2hsv(col2rgb(col.control$focus)))
            if(nbcol > 0) {
                bcol <- matrix(c(rep(0, 2*nbcol), rep(.9, nbcol)), ncol = 3)
                rbind(foc.col, bcol)
            }
            else foc.col
        }
    colnames(fills) <- c("h","s","v")
    diff.cols[[1]] <- list(fill = fills, border = gray(.8))

    ##_______________ Full Cell Plotting for Unique Bin1 Cells_________________

    if(length(levcells) > 0) {
        for(i in unique(whichbin)) {
            pcells <-
                if(eroded[i])
                    bin1@hbins[[i]]@cell[bin1@hbins[[i]]@eroded]
                else bin1@hbins[[i]]@cell
            pcells <- which(pcells %in% levcells[whichbin == i])
            pfill <- diff.cols[[1]]$fill[i,]
            pfill <- hsv(pfill[1],pfill[2],pfill[3])
            hexpolygon(x = tmph.xy[[i]]$x[pcells],
                       y = tmph.xy[[i]]$y[pcells], hexC,
                       border = diff.cols[[1]]$border ,
                       fill = pfill)
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
                              function(x) paste(which(x), sep = "", collapse = ":"))
            inter.nm <- unique(whichbin)
            #fills <- matrix(0, length(inter.nm), 3)
            fills <- rep(hsv(1), length(inter.nm))
            i <- 1
            for(bn in inter.nm) {
                who <- as.integer(unlist(strsplit(bn, ":")))
                fills[i] <- mixcolors2(diff.cols[[1]]$fill[who,],
                                       1/length(who),where = "LUV")
                i <- i+1
            }
            #fills <- LUV(fills)
            diff.cols[[dl]] <- list(fill = fills,
                                    border = gray((i.depth-dl)/i.depth))
            ##____Full Cell Plotting for Intersecting Cells at Intersection Depth i____
            i <- 1
            for(ints in inter.nm) {
                bin.i <- as.integer(unlist(strsplit(ints, ":"))[1])
                pcells <-
                    if(eroded[bin.i])
                        bin1@hbins[[bin.i]]@cell[bin1@hbins[[bin.i]]@eroded]
                    else bin1@hbins[[bin.i]]@cell
                pcells <- which(pcells %in% levcells[whichbin == ints])
                hexpolygon(x = tmph.xy[[bin.i]]$x[pcells],
                           y = tmph.xy[[bin.i]]$y[pcells], hexC,
                           border = diff.cols[[dl]]$border ,
                           fill = diff.cols[[dl]]$fill[i] )
                i <- i+1
            }
        }

    }

    ##_____________________________Plot Median Cells___________________________

    ## With all these colors floating around I think it would be worth
    ## porting the 3d hexagon stuff to grid. Then it would be easier
    ## to distinguish the medians because they would stand out like
    ## little volcanoes :) NL
    if(any(eroded)) {
        hmeds <- matrix(unlist(lapply(bin1@hbins[eroded],
                                      function(x) unlist(getHMedian(x)))),
                        ncol = 2, byrow = TRUE)
        hexpolygon(x = hmeds[, 1], y = hmeds[, 2], hexC,
                   border = col.control$med.b, fill = col.control$medhex)
        if(arrows) {
            for(i in focus) {
                for(j in focus[focus < i]) {
                    if(abs(hmeds[i, 1] - hmeds[j, 1]) +
                       abs(hmeds[i, 2] - hmeds[j, 2]) > eps)
			grid.lines(c(hmeds[i, 1],hmeds[j, 1]),
                                   c(hmeds[i, 2], hmeds[j, 2]),
                                   default.units = "native",
                                   arrow=arrow(length=size))
                        #grid.arrows(c(hmeds[i, 1], hmeds[j, 1]),
                        #            c(hmeds[i, 2], hmeds[j, 2]),
                        #            default.units = "native",
                        #            length = size, gp = gpar(lwd = lwd))
                }
            }
        }
    }

    ##________________Clean Up_______________________________________________

    popViewport()
    invisible(hvp)
} ## hdiffplot()
