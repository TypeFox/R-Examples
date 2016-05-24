plot.fdt_cat.default <- function (x,
                                  type=c('fb', 'fp', 'fd',
                                         'rfb', 'rfp', 'rfd',
                                         'rfpb', 'rfpp', 'rfpd',
                                         'cfb', 'cfp', 'cfd',
                                         'cfpb', 'cfpp', 'cfpd',
                                         'pa'),
                                  v=FALSE,
                                  v.round=2,
                                  v.pos=3,
                                  xlab=NULL,
                                  xlas=0,
                                  ylab=NULL,
                                  y2lab=NULL,
                                  y2cfp=seq(0, 100, 25),
                                  col=gray(.4),
                                  xlim=NULL,
                                  ylim=NULL,
                                  main=NULL,
                                  box=FALSE, ...)
{
  # barplot
  plot_b <- function(...)
  {
    cw <- max(sapply(names(y),
                     nchar))

    if (xlas == 1) 
      mar <- c(1, 
               1, 
               0, 
               2)
    else 
      mar <- c(log(max(cw),
                   2), 
               0, 
               0, 
               2)

    oldpar <- par(mar=pmax(par('mar') + mar, 
                           c(4.1, 
                             4.1, 
                             3.1, 
                             4.1)), 
                  las=xlas, 
                  no.readonly=TRUE)

    on.exit(par(oldpar))

    mp <-  barplot(y,
                   xlab=xlab,
                   col=col,
                   xlim=xlim,
                   ylim=ylim,   
                   main=main, 
                   yaxt='n', ...)

    rect(mp - 0.5,
         rep(0, 
             length(y)), 
         mp + 0.5,
         y, 
         col=col)

    if (box)
      box()

    axis(2, 
         las=3)

    mtext(ylab, 
          2, 
          line=2.5, 
          las=3)

    if (v) {
      text(mp,
           y,
           labels=format(round(y,
                               v.round),
                         nsmall=v.round),
           pos=v.pos,
           xpd=TRUE, 
           col=col, ...)
    }
  }

  # polygon
  plot_p <- function(...)
  {
    plot(1:dim(x)[1],
         y,
         type='b',
         xaxt='n',
         bty='n',
         xlim=xlim,
         ylim=ylim,
         xlab=xlab,
         ylab=ylab,
         col=col,
         main=main, ...)

    axis(side=1,
         at=1:dim(x)[1],
         labels=x$Category,
         las=xlas, ...)

    if (v) {
      text(1:dim(x)[1],  
           y,
           labels=format(round(y,
                               v.round),
                         nsmall=v.round),
           pos=v.pos,
           xpd=TRUE, 
           col=col, ...)
    }
  }

  # dotchart
  plot_d <- function(...)
  {
    if (is.null(xlim)) {
      ix <- min(y) * .1
      xmin <- min(y) - ix
      xmax <- max(y) + ix
      xlim <- c(xmin, xmax)
    }

    dotchart(y, 
             xlab=xlab,
             ylab=ylab,
             color=col,
             xlim=xlim,
             main=main, ...)

    if (v) {
      text(y,
           1:dim(x)[1],  
           labels=format(round(y,
                               v.round),
                         nsmall=v.round),
           pos=v.pos,
           xpd=TRUE, 
           col=col, ...)
    }
  }

  # pareto
  plot_pa <- function(...)
  {
    cw <- max(sapply(names(y),
                     nchar))

    if (xlas == 1) 
      mar <- c(1, 
               1, 
               0, 
               2)
    else 
      mar <- c(log(max(cw),
                   2), 
               0, 
               0, 
               2)

    oldpar <- par(mar=pmax(par('mar') + mar, 
                           c(4.1, 
                             4.1, 
                             3.1, 
                             4.1)), 
                  las=xlas, 
                  no.readonly=TRUE)

    on.exit(par(oldpar))

    mp <-  barplot(y,
                   xlab=xlab,
                   col=col,
                   xlim=xlim,
                   ylim=ylim,   
                   main=main, 
                   yaxt='n', ...)

    if (v) {
      # f
      vf <- c(paste(y[1], 
                    round(cfp[1],
                          1),
                    sep='-'),
              y[2:length(y)])

      text(mp,
           y,
           vf,
           pos=v.pos,
           xpd=TRUE, 
           col=col, ...)

      # cf-cfp
      vcfp <- c(paste(round(cf[-1],
                            1),
                      round(cfp[-1],
                            1),
                      sep='-'))

      ifelse(length(col) > 1,
             col.tmp <- col[-1],
             col.tmp <- col)

      text(mp[-1],
           cf[-1],
           vcfp,
           pos=v.pos,
           xpd=TRUE, 
           col=col.tmp, ...)
    }

    q <- quantile(seq(0, 
                      max(cf),
                      by=max(cf) / 100),
                  y2cfp / 100)

    abline(h=q, 
           col=gray(.8), 
           lty=3, ...)

    # clear lines over bars
    rect(mp - 0.5,
         rep(0, 
             length(y)), 
         mp + 0.5,
         y, 
         col=col, ...)

    lines(mp,
          cf, 
          type='b', 
          cex=0.7, 
          pch=19,
          col=col, ...)

    if (box)
      box()

    axis(2, 
         las=3, ...)

    mtext(ylab, 
          2, 
          line=2.5, 
          las=3, ...)

    axis(4, 
         at=q, 
         las=3, 
         labels=round(y2cfp,
                      2), ...)

    mtext(y2lab,
          4, 
          line=2.5, 
          las=3, ...)
  }

  switch(match.arg(type),
         # fb (frequency) - bar
         fb = {
           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- x[, 2]
           names(y) <- x$Category

           if (is.null(ylim))
             ylim=c(0, 
                    max(y) * 1.3)

           plot_b(...)
         },

         # fp (frequency) - polygon
         fp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(x[, 2])))

           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- x[, 2]

           plot_p(...)
         },

         # fd (frequency) - dotchart
         fd = {
           y <- x[, 2]
           names(y) <- x$Category

           if(is.null(xlab))
             xlab <- 'Frequency'

           plot_d(...)
         },

         # pa (frequency) - pareto
         pa = {
           y <- x[, 2]
           names(y) <- x$Category

           cf <- x[, 5]
           cfp <- x[, 6]

           if (is.null(xlab))
             xlab <- 'Category'
           if (is.null(ylab))
             ylab <- 'Frequency'
           if (is.null(y2lab))
             y2lab <- 'Cumulative frequency, (%)'

           if (is.null(ylim))
             ylim=c(0, 
                    sum(y) * 1.1)

           plot_pa(...)
         },

         # rfb (relative frequency) - bar
         rfb = {
           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Relative frequency'

           y <- x[, 3]
           names(y) <- x$Category

           plot_b(...)
         },  

         # rfp (relative frequency) - polygon
         rfp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(x[, 3])))

           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Relative frequency'

           y <- x[, 3]

           plot_p(...)
         },

         # rfd (relative frequency) - dotchart
         rfd = {
           y <- x[, 3]
           names(y) <- x$Category

           if(is.null(xlab))
             xlab <- 'Relative frequency'

           plot_d(...)
         },

         # rfpb (relative frequency %) - bar
         rfpb = {
           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Relative frequency (%)'

           y <- x[, 4]
           names(y) <- x$Category

           plot_b(...)
         },  

         # rfpp (relative frequency %) - polygon
         rfpp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(x[, 4])))

           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Relative frequency (%)'

           y <- x[, 4]

           plot_p(...)
         },

         # rfpd (relative frequency %) - dotchart
         rfpd = {
           y <- x[, 4]
           names(y) <- x$Category

           if(is.null(xlab))
             xlab <- 'Relative frequency (%)'

           plot_d(...)
         },

         # cfb (cumulative frequency) - bar
         cfb = {
           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Cumulative frequency'

           y <- x[, 5]
           names(y) <- x$Category

           plot_b(...)
         },  

         # cfp (cumulative frequency) - polygon
         cfp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(x[, 5])))

           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Cumulative frequency'

           y <- x[, 5]

           plot_p(...)
         },

         # cfd (cumulative frequency) - dotchart
         cfd = {
           y <- x[, 5]
           names(y) <- x$Category

           if(is.null(xlab))
             xlab <- 'Cumulative frequency'

           plot_d(...)
         },

         # cfpb (cumulative frequency %) - bar
         cfpb = {
           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Cumulative frequency (%)'

           y <- x[, 6]
           names(y) <- x$Category

           plot_b(...)
         },  

         # cfpp (cumulative frequency %) - polygon
         cfpp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(x[, 6])))

           if(is.null(xlab))
             xlab <- 'Category'
           if(is.null(ylab))
             ylab <- 'Cumulative frequency (%)'

           y <- x[, 6]

           plot_p(...)
         },

         # cfpd (cumulative frequency %) - dotchart
         cfpd = {
           y <- x[, 6]
           names(y) <- x$Category

           if(is.null(xlab))
             xlab <- 'Cumulative frequency (%)'

           plot_d(...)
         })
}
