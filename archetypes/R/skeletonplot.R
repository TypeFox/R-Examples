

#' Skeleton plot.
#'
#' Displays a schematic representation of skeleton data as available
#' in dataset \code{\link{skel}}.
#'
#' @param x Matrix or data.frame of skeleton data.
#' @param skel.width Reference width for instance calculation.
#' @param skel.height Reference height for instance calculation.
#' @param base.radius Base radius for points.
#' @param xlab The x label of the plot.
#' @param ylab The y label of the plot.
#' @param xlim Numeric of length 2 giving the x limits for the plot.
#' @param ylim Numeric of length 2 giving the y limits for the plot.
#' @param col Color of the different parts of the skeleton.
#' @param mtext Label archetypes.
#' @param skel.lwd Line width of skeleton.
#' @param ... Passed to underlying canvas plot function.
#' @return List of skeleton instances.
#' @export
#' @seealso \code{\link{skel}}
skeletonplot <- function(x, skel.width = 100, skel.height = 200,
                         ylab = 'Height (cm)', base.radius = 2, xlab = '',
                         xlim = (nrow(x)*c(0,skel.width)), ylim = c(0, skel.height),
                         col = NULL, mtext = TRUE, skel.lwd = 1, ...) {

  if ( is.data.frame(x) )
    x <- as.matrix(x)
  
  if ( is.null(col) ) {
    col <- c(hipbase = 1, hip = 1, shoulderbase = 1, shoulder = 1,
      head = 1, elbow = 2, wrist = 3, knee = 4, ankle = 5,
      chest = 'purple1', pelvis = 6)
  }

  ### Skeleton model (see human-modelling.vsd):
  model.y <- c(ankle=0, knee=7, wrist=12, hip=13, hipbase=15, pelvis=16,
               waist=17, elbow=20, chest=24, shoulder=26,
               shoulderbase=27, head=30, top=32) / 32

  model.x.leg <- c(hip=1, knee=1.5, ankle=1)
  model.x.spine <- c(hipbase=0, pelvis=0, waist=0, chest=0, shoulderbase=0, head=0, top=0)
  model.x.arm <- c(shoulder=1, elbow=5/3, wrist=5/3)


  ### One skeleton instance:
  one.skeleton <- function(x, x0=0) {

    # Calculate instance:
    skel.y <- model.y * x['Height']

    skel.x.leg <- model.x.leg * (x['Bitro'] / 2)
    skel.x.spine <- model.x.spine
    skel.x.arm <- model.x.arm * (x['Biac'] / 2)
    skel.x <- c(skel.x.leg, skel.x.spine, skel.x.arm)

    skel.circles <- base.radius + c(hipbase=0, hip=0, shoulderbase=0,
                                    shoulder=0, head=0,
                                    elbow=unname(x['ElbowDiam'])/2,
                                    wrist=unname(x['WristDiam'])/2,
                                    knee=unname(x['KneeDiam'])/2,
                                    ankle=unname(x['AnkleDiam'])/2) / 2

    skel.rectangles <- rbind(chest=c(width=unname(x['ChestDiam']),
                               height=unname(x['ChestDp'])),
                             pelvis=c(width=unname(x['Biil']), height=0))


    # Plot it:
    lines(x0 + skel.x.spine, skel.y[names(skel.x.spine)],
          lwd=skel.lwd, ...)

    lines(x0 + c(skel.x.spine['hipbase'], skel.x.leg),
          c(skel.y['hipbase'], skel.y[names(skel.x.leg)]),
          lwd=skel.lwd, ...)
    lines(x0 + c(skel.x.spine['hipbase'], -skel.x.leg),
          c(skel.y['hipbase'], skel.y[names(skel.x.leg)]),
          lwd=skel.lwd, ...)

    lines(x0 + c(skel.x.spine['shoulderbase'], skel.x.arm),
          c(skel.y['shoulderbase'], skel.y[names(skel.x.arm)]),
          lwd=skel.lwd, ...)
    lines(x0 + c(skel.x.spine['shoulderbase'], -skel.x.arm),
          c(skel.y['shoulderbase'], skel.y[names(skel.x.arm)]),
          lwd=skel.lwd, ...)

    d <- names(skel.circles)
    symbols(x0 + c(skel.x[d], -skel.x[d]), rep(skel.y[d], 2),
            circles=rep(skel.circles[d], 2), fg=col[d], bg=col[d],
            xlim=xlim, ylim=ylim, inches=FALSE, add=TRUE)

    d <- rownames(skel.rectangles)
    symbols(rep(x0, length(d)), skel.y[d],
            rectangles=skel.rectangles, fg=col[d],
            xlim=xlim, ylim=ylim, inches=FALSE, add=TRUE, lwd=2)

    return(list(x0=x0, x=skel.x, y=skel.y,
                circles=skel.circles, rectangles=skel.rectangles))
  }


  ### Plot:
  nskels <- nrow(x)

  yticks <- seq(ylim[1], ylim[2], by=20)
  xticks <- seq(xlim[1], xlim[2], by=50)

  # Canvas:
  plot(1, xlim=xlim, ylim=ylim, type='n', xlab=xlab, ylab=ylab, axes=FALSE, ...)
  axis(1, at=xticks)
  axis(2, at=yticks)
  box()

  # Gridlines:
  abline(v=xticks, col='lightgray', lty='dotted', lwd=1)
  abline(h=yticks, col='lightgray', lty='dotted', lwd=1)

  # Skeletons:
  skels <- list()

  for ( i in 1:nskels ) {
    x0 <- (i-1) * skel.width + (skel.width/2)

    skels[[i]] <- one.skeleton(x[i,], x0=x0)

    if ( mtext )
      mtext(paste('Archetype', i), side=3, line=0, at=x0)
  }


  invisible(skels)
}



#' Annotated skeleton plot.
#'
#' Displays a generic skeleton with annotations explaining the
#' measurements available in data set \code{\link{skel}}.
#'
#' @return Generic skeleton instance.
#' @rdname skeletonplot
#'
#' @export
jd <- function() {
  jd <- rbind(c(AnkleDiam=13.9, KneeDiam=18.8, WristDiam=10.5, Bitro=32.0,
                Biil=27.8, ElbowDiam=13.4, ChestDiam=28.0, ChestDp=15,
                Biac=38.8, Height=171.1))

  s <- skeletonplot(jd, skel.height=190,
                    mtext=FALSE, xlim=c(-100,200), skel.lwd=1)[[1]]


  ### Annotate JD:
  acol <- gray(0.5)

  annotation1 <- function(text, x, y, alen=10) {
    ws <- 0

    arrows(x, y, x+alen, y,
           length=0.1, code=1, col=acol, lwd=1)

    text(labels=text,
         x=x+alen+ws, y=y, pos=ifelse(alen<0,2,4), col=acol)
  }

  annotation2 <- function(text, xb, xd, y, offset, alen=-30) {
    x0 <- xb - xd
    x1 <- xb + xd

    lines(c(x0, x0), c(y, y+offset), col=acol)
    lines(c(x1, x1), c(y, y+offset), col=acol)

    arrows(x0, y+offset, x1, y+offset,
           code=3, length=0.1, col=acol)

    lines(c(x0, x0+alen), c(y+offset, y+offset), col=acol)

    text(labels=text, x=x0+alen, y=y+offset, pos=2, col=acol)
  }


  annotation1('Diameter of ankle',
              s$x0 + s$x['ankle'] + s$circles['ankle'],
              s$y['ankle'])

  annotation1('Diameter of knee',
              s$x0 + s$x['knee'] + s$circles['knee'],
              s$y['knee'])

  annotation1('Diameter of wrist',
              s$x0 + s$x['wrist'] + s$circles['wrist'],
              s$y['wrist'])

  annotation1('Diameter of pelvis\n(biiliac)',
              s$x0 + s$x['pelvis'] + s$rectangles['pelvis','width']/2,
              s$y['pelvis'], alen=25)

  annotation1('Diameter of elbow',
              s$x0 + s$x['elbow'] + s$circles['elbow'],
              s$y['elbow'])

  annotation1('Height',
              s$x0 + s$x['top'],
              s$y['top'])

  annotation2('Diameter between\nhips (bitrochanteric)',
              s$x0, s$x['hip'], s$y['hip'] - s$circles['hip'], -15)

  annotation2('Diameter between\nshoulders (biacromial)',
              s$x0, s$x['shoulder'], s$y['shoulder'] + s$circles['shoulder'], 10)

  annotation2('Diameter of chest',
              s$x0, s$rectangles['chest','width']/2,
              s$y['chest'] - s$rectangles['chest','height']/2, -5)

  annotation1('Depth of chest',
              s$x0 - s$x['chest'] - s$rectangles['chest','width']/2,
              s$y['chest'], alen=-20)

  invisible(s)
}

