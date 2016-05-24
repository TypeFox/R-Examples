## Plot objects of the class 'bpca.3d' with the packages 'scatterplot3d'
## ('graphics' based) and 'rgl'
plot.bpca.3d <- function(x,
                         rgl.use=FALSE,
                         ref.lines=TRUE,
                         ref.color='navy',
                         ref.lty=ifelse(rgl.use, NA, 'dotted'),
                         clear3d=ifelse(rgl.use, TRUE, NULL),
                         simple.axes=ifelse(rgl.use, TRUE, NULL),
                         aspect=ifelse(rgl.use, c(1, 1, 1), NULL),
                         var.factor=1,
                         var.color='red3',
                         var.lty=ifelse(rgl.use, NA, 'solid'), 
                         var.pch=ifelse(rgl.use, NULL, 20),
                         var.pos=ifelse(rgl.use, 0, 4),
                         var.cex=ifelse(rgl.use, .8, .6),
                         var.offset=ifelse(rgl.use, NULL, .2),
                         obj.color='black',
                         obj.pch=ifelse(rgl.use, NULL, 20),
                         obj.pos=ifelse(rgl.use, 0, 4),
                         obj.cex=ifelse(rgl.use, .8, .6),
                         obj.offset=ifelse(rgl.use, NULL, .2),
                         obj.names=TRUE,
                         obj.labels,
                         obj.identify=FALSE,
                         box=FALSE,
                         angle=ifelse(rgl.use, NULL, 40),
                         xlim, ylim, zlim, xlab, ylab, zlab, ...)
{
  if (!inherits(x, 
                'bpca.3d'))
    stop("Use this function only with 'bpca.3d' class!")

  d1 <- x$number[1]; d3 <- x$number[3]
  coobj <- x$coord$objects[,d1:d3]
  covar <- x$coord$variables[,d1:d3]

  scores <- rbind(coobj,
                  covar * var.factor,
                  rep(0, 
                      ncol(coobj)))

  # spatial limits
  if (missing(xlim) || missing(ylim) || missing(zlim)) {
    ms  <- max(abs(scores)) * 1.2
    msp <- c(-ms, ms)
  }  

  if (missing(obj.labels))
    obj.labels <- rownames(coobj)

  if (missing(xlim))
    xlim <- msp
  if (missing(ylim))
    ylim <- msp
  if (missing(zlim))
    zlim <- msp

  if (missing(xlab) || missing(ylab) || missing(zlab)) {
    eigv <- x$eigenvalues
    prop <- 100 * eigv^2 / sum(eigv^2)
    labs <- paste('PC',
                  d1:d3,
                  ' (',
                  round(prop[d1:d3], 
                        2),
                  '%)', 
                  sep='')
  }

  if (missing(xlab))
    xlab <- labs[1]
  if (missing(ylab))
    ylab <- labs[2]
  if (missing(zlab))
    zlab <- labs[3]

  # Plot bpca.3d under package 'scatterplot3d'
  if(!rgl.use) {
    op <- par(no.readonly=TRUE)
    # a empty plot (reference)
    graph <- scatterplot3d(scores,
                           xlim=xlim,
                           ylim=ylim,
                           zlim=zlim,
                           type='n',
                           xlab=xlab,
                           ylab=ylab,
                           zlab=zlab,
                           grid=FALSE,
                           box=box,
                           angle=angle, ...)

    # objects
    if(obj.names) {
      # points of objects
      graph$points3d(coobj,
                     pch=obj.pch,
                     type='p',
                     col=obj.color,
                     cex=obj.cex, ...)

      # labels of objects
      text(graph$xyz.convert(coobj),
           labels=obj.labels,
           pos=obj.pos,
           offset=obj.offset,
           col=obj.color,
           cex=obj.cex, ...)
    } else {
      graph$points3d(coobj,
                     pch=obj.pch,
                     type='p',
                     col=obj.color,
                     cex=obj.cex, ...)
    }

    # variables
    for(i in 1:nrow(covar)) {
      # points of variables
      graph$points3d(c(0, 
                       covar[i,1] * var.factor),
                     c(0, 
                       covar[i,2] * var.factor),
                     c(0, 
                       covar[i,3] * var.factor),
                     pch=var.pch,
                     col=var.color,
                     type='p',
                     lty=var.lty,
                     cex=var.cex, ...)

      # segments of variables (vectors)
      graph$points3d(c(0, 
                       covar[i,1] * var.factor),
                     c(0, 
                       covar[i,2] * var.factor),
                     c(0, 
                       covar[i,3] * var.factor),
                     col=var.color,
                     type='l',
                     lty=var.lty, ...)
    }

    # labels of variables
    text(graph$xyz.convert(covar * var.factor),
         labels=rownames(covar),
         pos=var.pos,
         offset=var.offset,
         col=var.color,
         cex=var.cex, ...)

    # reference lines
    if(ref.lines) {
      # x
      graph$points3d(xlim,
                     c(0, 0),
                     c(0, 0),
                     type='l',
                     lty=ref.lty,
                     col=ref.color, ...)

      # y
      graph$points3d(c(0, 0),
                     ylim,
                     c(0, 0),
                     type='l',
                     lty=ref.lty,
                     col=ref.color, ...)

      # z
      graph$points3d(c(0, 0),
                     c(0, 0),
                     zlim,
                     type='l',
                     lty=ref.lty,
                     col=ref.color, ...)
    }

    if(obj.identify) 
      identify(x=graph$xyz.convert(coobj),
               labels=obj.labels,
               cex=obj.cex)

    par(op,
        no.readonly=TRUE)
  }

  # Plot bpca.3d under package 'rgl'
  if(rgl.use) {
    size <- max(coobj) /
    20*obj.cex

    if (clear3d)
      clear3d()

    # a empty plot (reference)
    plot3d(scores,
           xlim=xlim,
           ylim=ylim,
           zlim=zlim,
           xlab='',
           ylab='',
           zlab='',
           type='n',
           axes=FALSE,
           box=box,
           aspect=aspect,
           top=TRUE, ...)

    # objects
    if (obj.names) {
      # points of objects
      spheres3d(coobj,
                col=obj.color,
                radius=size / 2,
                alpha=.5, ...)

      # labels of objects
      text3d(coobj,
             texts=obj.labels,
             col=obj.color,
             alpha=.5,
             adj=obj.pos,
             cex=obj.cex, ...)
    } else {
      spheres3d(coobj,
                col=obj.color,
                radius=size,
                alpha=.5, ...)
    }
    # variables
    for(i in 1:nrow(covar)) {
      # points of variables
      spheres3d(covar[i,] * var.factor,
                col=var.color,
                radius=size / 2,
                alpha=.5, ...)

      # segments of variables (vectors)
      segments3d(rbind(matrix(0,
                              ncol=3),
                       covar[i,] * var.factor),
                 col=var.color, ...)
    }

    # labels of variables
    text3d(covar * var.factor,
           texts=rownames(covar),
           col=var.color,
           adj=var.pos,
           cex=var.cex, ...)

    # axes
    if(simple.axes) {
      axes3d(c('x', 'y', 'z'))
      # simple axis
      title3d(xlab=xlab,
              ylab=ylab,
              zlab=zlab, ...)
    } else
      # axes with box
      decorate3d(xlab=xlab,
                 ylab=ylab,
                 zlab=zlab,
                 box=box, ...)

    # reference lines
    if(ref.lines) {
      # x
      lines3d(xlim,
              c(0, 0),
              c(0, 0),
              lty=ref.lty,
              col=ref.color, ...)

      # y
      lines3d(c(0, 0),
              ylim,
              c(0, 0),
              lty=ref.lty,
              col=ref.color, ...)

      # z
      lines3d(c(0, 0),
              c(0, 0),
              zlim,
              lty=ref.lty,
              col=ref.color, ...)
    }
  }
}

