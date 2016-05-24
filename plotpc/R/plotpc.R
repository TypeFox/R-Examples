# plotpc.R
# TODO position of x axis and y axis histograms is sometimes wrong

plotpc <- function(x, # x is a two column matrix or data.frame
    xrange=NULL,
    hist=TRUE,
    main="Principal components",
    xlab=NULL,
    ylab=NULL,
    gp.points=gpar(cex=.6),
    pch=20,
    height=xrange/10,
    breaks="Sturges", # "breaks" param for hist (default same as hist function)
    adjust=1,         # "adjust" param for density
    gp.hist=if(hist) gp.hist <- gpar(col="gray", fill="gray")
            else     gp.hist <- gpar(col="black"),
    gp.text=gpar(cex=.8, font=2),
    gp.axis=gpar(col="gray", lwd=2),
    sd.ellipse=NA,
    gp.ellipse=gpar(col="gray", lwd=2),
    heightx=NULL, breaksx=NULL, adjustx=NULL, gp.histx=NULL,
                 textx="", gp.textx=NULL, axis.lenx=0, gp.axisx=NULL,
    heighty=NULL, breaksy=NULL, adjusty=NULL, gp.histy=NULL,
                 texty="", gp.texty=NULL, axis.leny=0, gp.axisy=NULL,
    height1=NULL, flip1=FALSE,
                 breaks1=NULL, adjust1=NULL, gp.hist1=NULL, offset1=NULL,
                 text1=NULL, gp.text1=NULL, axis.len1=2, gp.axis1=NULL,
    height2=NULL, flip2=FALSE,
                 breaks2=NULL, adjust2=NULL, gp.hist2=NULL, offset2=NULL,
                 text2=NULL, gp.text2=NULL, axis.len2=2, gp.axis2=NULL,
    angle3=NA, height3=NULL, flip3=FALSE,
                 breaks3=NULL, adjust3=NULL, gp.hist3=NULL, offset3=NULL,
                 text3=NULL, gp.text3=NULL, axis.len3=0, gp.axis3=NULL,
    angle4=NA, height4=NULL, flip4=FALSE,
                 breaks4=NULL, adjust4=NULL, gp.hist4=NULL, offset4=NULL,
                 text4=NULL, gp.text4=NULL, axis.len4=0, gp.axis4=NULL,
    angle5=NA, height5=NULL, flip5=FALSE,
                 breaks5=NULL, adjust5=NULL, gp.hist5=NULL, offset5=NULL,
                 text5=NULL, gp.text5=NULL, axis.len5=0, gp.axis5=NULL,
    angle6=NA, height6=NULL, flip6=FALSE,
                 breaks6=NULL, adjust6=NULL, gp.hist6=NULL, offset6=NULL,
                 text6=NULL, gp.text6=NULL, axis.len6=0, gp.axis6=NULL,
    angle7=NA, height7=NULL, flip7=FALSE,
                 breaks7=NULL, adjust7=NULL, gp.hist7=NULL, offset7=NULL,
                 text7=NULL, gp.text7=NULL, axis.len7=0, gp.axis7=NULL,
    yonx = FALSE, offset.yonx=-xrange/2.5,
                 text.yonx="y~x", gp.text.yonx=NULL,
                 axis.len.yonx=xrange/2.5, gp.axis.yonx=gpar(col=1),
    xony = FALSE, offset.xony=-xrange/2.5,
                 text.xony="x~y", gp.text.xony=NULL,
                 axis.len.xony=xrange/2.5, gp.axis.xony=gpar(col=1))
{
    to.radians <- function(angle) angle <- angle * pi / 180

    to.degrees <- function(angle) angle <- 180 * angle / pi

    is.valid <- function(angle)
         !is.na(angle) && !is.null(angle) && is.numeric(angle) && length(angle) == 1

    vp1 <- function(xscale, yscale, angle=0)
        viewport(width=unit(.7, "npc"), height=unit(.7, "npc"), # .7 for whitespace
                 xscale=xscale, yscale=yscale, angle=angle)

    # get angle to rotate the text by, either 0 or 180, depending on angle1
    # i.e. draw text the right way up or upside down, relative to the viewport
    # note: pass in radians but returns degrees
    get.rot <- function(angle1) {
        if(angle1 < 0)
            angle1 <- angle1 + 2 * pi
        if(angle1 >= 0 && angle1 <= pi/2)
            180
        else if(angle1 > 3/2 * pi && angle1 < 2 * pi)
            180
        else
            0
    }
    plot1 <- function(angle1, flip1, height1, breaks1, adjust1, gp.hist1, offset1,
                      axis.len1, gp.axis1, text1, text1.default, gp.text1) {
        # get defaults from globals: x, xlim, height, breaks, adjust, ...
        if(is.null(height1))
            height1 <- height
        if(is.null(breaks1))
            breaks1 <- breaks
        if(is.null(adjust1))
            adjust1 <- adjust
        if(is.null(gp.hist1))
            gp.hist1 <- gp.hist
        if(is.null(offset1))
            offset1 <- offset.default
        if(is.null(gp.axis1))
            gp.axis1 <- gp.axis
        if(is.null(text1))
            text1 <- text1.default
        if(is.null(gp.text1))
            gp.text1 <- gp.text
        if(angle1 >= 2 * pi)
            angle1 <- angle1 - 2 * pi
        if(angle1 <= -2 * pi)
            angle1 <- angle1 + 2 * pi
        angle1 <- angle1 + pi   # TODO why?
        x.rot <- x[,1] * cos(angle1) + x[,2] * sin(angle1) # rotated observations
        pushViewport(vp1(xlim, xlim, to.degrees(angle1)))
        if(axis.len1 > 0) {     # draw axis?
            len1 <- sd(x.rot) * axis.len1
            grid.lines(x=un(c(-len1, len1)), y=un(c(0,0)), gp=gp.axis1)
        }
        if(height1) {
            if(hist)
                my.plot.histogram(hist(x.rot, breaks=breaks1, plot=FALSE),
                                  height1, flip1, gp.hist1, offset1)
            else
                my.plot.density(density(x.rot, adjust=adjust1),
                                height1, flip1, gp.hist1, offset1)

            if(nchar(text1)) {  # draw text1 above histogram?
                y1 <- offset1 + max(0, height1) + (xlim[2] - xlim[1]) / 20
                if(flip1)
                    y1 <- -y1
                # TODO would better spacing with "lines" not "native" for y1 below?
                grid.text(text1, x=un(0), y=un(y1), gp=gp.text1,
                          rot=get.rot(angle1-pi))
            }
        }
        popViewport(1)
    }
    plot.yonx <- function(angle, offset,
                          axis.len, gp.axis, text, text.default, gp, above) {
        if(is.null(gp))
            gp <- gp.text
        pushViewport(vp1(xlim, xlim, to.degrees(angle)))
        # draw the axis
        grid.lines(x=un(c(-axis.len, axis.len)), y=un(c(0,0)), gp=gp.axis)
        if(nchar(text))
            grid.text(text, just="left", x=un(offset),
                      y=un(if(above) xrange/30 else -xrange/50),
                      gp=gp, rot=get.rot(angle + pi))
        popViewport(1)
    }
    # plotpc starts here

    if(missing(x))
        stop("x argument is missing")
    x <- as.matrix(x)
    if(NCOL(x) < 2)
        stop("x must have 2 columns")
    if(NCOL(x) > 2)
        warning("x has more than 2 columns, using only columns 1 and 2")
    if(NROW(x) < 3) # 3 is arbitrary
        stop("x has too few rows")
    if(is.null(colnames(x)))
        colnames(x) <- c("x1", "x2")
    if(is.null(xlab))
        xlab <- colnames(x)[1]
    if(is.null(ylab))
        ylab <- colnames(x)[2]

    xmean <- mean(x[,1], na.rm=FALSE)
    if(is.na(xmean))
        stop("NA in x[,1]")
    ymean <- mean(x[,2], na.rm=FALSE)
    if(is.na(ymean))
        stop("NA in x[,2]")
    x <- scale(x, scale=FALSE) # mean-center x

    if(is.null(xrange))
        xrange <- 4 * max(abs(x[,1:2]))
    xlim <- c(-xrange/2, xrange/2)
    if(xrange <= 0)
        stop("xrange <= 0")

    pc <- princomp(x[,1:2], scores=TRUE)
    offset.default <- xrange / 3
    grid.newpage()

    # plot confidence ellipse (we plot it first so it is under the points)

    if(is.valid(sd.ellipse) && sd.ellipse > 1) {
        angle <- atan2(pc$loadings[2,1], pc$loadings[1,1])
        pushViewport(vp1(xlim, xlim, to.degrees(angle)))
        a <- seq(0, 2 * pi, len=300)
        grid.lines(un(sd.ellipse * sd(pc$scores[,1]) * cos(a)),
                   un(sd.ellipse * sd(pc$scores[,2]) * sin(a)),
                   gp=gp.ellipse)
        popViewport(1)
    }
    # draw title, axes, and points (with the non-mean-centered x)

    vp <- vp1(xmean+xlim, ymean+xlim)
    pushViewport(vp)
    grid.text(main, y=un(ymean+xlim[2]+xrange/10), gp=gpar(cex=1.0, font=2))
    grid.xaxis(gp=gpar(cex=.8))
    grid.text(xlab, y=unit(-2.5, "lines"), gp=gpar(cex=.8))
    grid.yaxis(gp=gpar(cex=.8))
    grid.text(ylab, x=unit(-3.5, "lines"), gp=gpar(cex=.8), rot=90)
    grid.points(x[,1]+xmean, x[,2]+ymean, gp=gp.points, pch=pch)
    popViewport(1)

    # plot histogram for x[,1] on main x axis

    plot1(0, TRUE, heightx, breaksx, adjustx, gp.histx, xlim[1] + xrange / 100,
          axis.lenx, gp.axisx, textx, "", gp.textx)

    # plot histogram for x[,2] on main y axis

    plot1(pi/2, FALSE, heighty, breaksy, adjusty, gp.histy, -xlim[2] + xrange / 100,
          axis.leny, gp.axisy, texty, "", gp.texty)

    # plot histograms for angle3 to angle7

    if(is.valid(angle3))
        plot1(to.radians(angle3), flip3, height3, breaks3, adjust3, gp.hist3, offset3,
              axis.len3, gp.axis3, text3, sprintf(" %.0f deg", angle3), gp.text3)

    if(is.valid(angle4))
        plot1(to.radians(angle4), flip4, height4, breaks4, adjust4, gp.hist4, offset4,
              axis.len4, gp.axis4, text4, sprintf(" %.0f deg", angle4), gp.text4)

    if(is.valid(angle5))
        plot1(to.radians(angle5), flip5, height5, breaks5, adjust5, gp.hist5, offset5,
              axis.len5, gp.axis5, text5, sprintf(" %.0f deg", angle5), gp.text5)

    if(is.valid(angle6))
        plot1(to.radians(angle6), flip6, height6, breaks6, adjust6, gp.hist6, offset6,
              axis.len6, gp.axis6, text6, sprintf(" %.0f deg", angle6), gp.text6)

    if(is.valid(angle7))
        plot1(to.radians(angle7), flip7, height7, breaks7, adjust7, gp.hist7, offset7,
              axis.len7, gp.axis7, text7, sprintf(" %.0f deg", angle7), gp.text7)

    # plot histogram for first principal component

    angle1 <- atan2(pc$loadings[2,1], pc$loadings[1,1])
    plot1(angle1, flip1, height1, breaks1, adjust1, gp.hist1, offset1,
          axis.len1, gp.axis1, text1,
          sprintf("%.0f deg (PC1)", to.degrees(angle1)), gp.text1)

    # plot histogram for second principal component

    angle2 <- atan2(pc$loadings[2,2], pc$loadings[1,2]) + pi
    plot1(angle2, flip2, height2, breaks2, adjust2, gp.hist2, offset2,
          axis.len2, gp.axis2, text2,
          sprintf("%.0f deg (PC2)", to.degrees(angle2)), gp.text2)

    # plot regression lines

    if(yonx || xony) {
        lm.yonx <- lm(x[,2] ~ x[,1])
        angle.yonx <- atan(lm.yonx$coefficients[2])
        lm.xony <- lm(x[,1] ~ x[,2])
        angle.xony <- atan(1/lm.xony$coefficients[2])
        # "above" specifies to plot the text above or below the regression line
        above <- angle.yonx < angle.xony
        if(yonx)
            plot.yonx(angle.yonx, offset.yonx, axis.len.yonx, gp.axis.yonx,
                      text.yonx, NULL, gp.text.yonx, above)
        if(xony)
            plot.yonx(angle.xony, offset.xony, axis.len.xony, gp.axis.xony,
                      text.xony, NULL, gp.text.xony, !above)
    }
    invisible(vp)  # return vp so can easily add text later in "native" coords
}
un <- function(x) unit(x, "native")

my.plot.histogram <-  function(x, height, flip, gp, offset) # x is the histogram
{
    rect <- function(xleft, ybottom, xright, ytop, gp) # args same as as graphics:rect
    {
        n <- max(length(xleft), length(xright), length(ybottom), length(ytop))
        x <- rbind(rep.int(NA, n), xleft, xright, xright, xleft)[-1]
        y <- rbind(rep.int(NA, n), ybottom, ybottom, ytop, ytop)[-1]
        grid.polygon(un(x), un(y), gp=gp)
    }
    breaks <- x$breaks
    counts <- height * x$counts / max(x$counts)
    if(flip) {
        offset <- -offset
        counts <- -counts
    }
    rect(breaks[-length(breaks)], offset, breaks[-1], offset + counts, gp)
}
my.plot.density <-  function(x, height, flip, gp, offset) # x is the density
{
    y <- height * x$y / max(x$y)
    if(flip) {
        offset <- -offset
        y <- -y
    }
    # draw the density graph
    grid.lines(un(x$x), un(y + offset), gp=gp)
    # draw the base line
    grid.lines(x=un(c(min(x$x),max(x$x))), y=un(c(offset,offset)), gp=gp)
}
