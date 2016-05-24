plot.fdt.default <- function (x,
                              type=c('fh', 'fp', 
                                     'rfh', 'rfp', 'rfph', 'rfpp',
                                     'd', 'cdh', 'cdp', 
                                     'cfh', 'cfp', 'cfph', 'cfpp'),
                              v=FALSE,
                              v.round=2,
                              v.pos=3,
                              xlab='Class limits',
                              xlas=0,
                              ylab=NULL,
                              col='gray',
                              xlim=NULL,
                              ylim=NULL,
                              main=NULL,
                              x.round=2, ...)
{
  brk <- with(x,
              seq(breaks['start'],
                  breaks['end'],
                  breaks['h']))

  if (is.null(xlim))
    xlim <- with(x,
                 c(breaks['start'],
                   breaks['end']))

  mids <- 0.5 * (brk[-1] + 
                 brk[-length(brk)])

  switch(match.arg(type),
         # f (absolute frequency) - histogram
         fh = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 2])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 2]

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)
           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },

         # f (absolute frequency) - polygon
         fp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 2])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- with(x,
                     table[, 2])

           plot(mids,
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
           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },

         # rf (relative frequency) - histogram
         rfh = {
           h <- with(x, 
                     breaks[3])
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 3])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab, 
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 3]

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # rf (relative frequency) - polygon
         rfp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 3])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- with(x,
                     table[, 3])
           plot(mids,
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

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },

         # rf (relative frequency %) - histogram
         rfph = {
           h <- with(x,
                     breaks[3])
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 4])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 4]

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # rf (relative frequency %) - polygon
         rfpp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 4])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- with(x,
                     table[, 4])
           plot(mids,
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

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },

         # Density
         d = {
           h <- with(x, 
                     breaks[3])

           if (is.null(ylim))
             ylim <- with(x, 
                          c(0, 
                            1.2 * max(table[, 3] / h)))

           if(is.null(ylab))
             ylab <- 'Density'


           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 3] / h

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # cd (cumulative density) - histogram
         cdh = {
           h <- with(x, 
                     breaks[3])

           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 1.2))

           if(is.null(ylab))
             ylab <- 'Cumulative density'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- cumsum(x$table[, 3])

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # cm (cumulative density) - polygon
         cdp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 1.2))

           if(is.null(ylab))
             ylab <- 'Cumulative density'

           y <- c(0,
                  cumsum(x$table[, 3]))
           plot(brk,
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

           if(v)
             text(x=brk,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         }, 

         # cf (cumulative frequency) - histogram
         cfh = {
           h <- with(x,
                     breaks[3])

           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * max(table[, 5])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 5]

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # cf (cumulative frequency) - polygon
         cfp = {
           if (is.null(ylim))
             ylim <- with(x,
                          c(0, 
                            1.2 * sum(table['f'])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- with(x,
                     c(0,
                       table[, 5]))
           plot(brk,
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

           if(v)
             text(x=brk,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },

         # cfp (cumulative frequency %) - histogram
         cfph = {
           h <- with(x,
                     breaks[3])

           if (is.null(ylim))
             ylim <- with(x, 
                          c(0, 
                            1.2 * max(table[, 6])))

           if(is.null(ylab))
             ylab <- 'Frequency'

           plot.new()

           plot.window(xlim,
                       ylim)

           title(main=main,
                 xlab=xlab,
                 ylab=ylab, ...)

           axis(2, ...)

           y <- x$table[, 6]

           rect(brk[-length(brk)],
                0,
                brk[-1],
                y,
                col=col, ...)

           if(v)
             text(x=mids,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         },  

         # cfp (cumulative frequency %) - polygon
         cfpp = {
           if (is.null(ylim))
             ylim <- c(0, 
                       1.2 * 100)

           if(is.null(ylab))
             ylab <- 'Frequency'

           y <- with(x,
                     c(0, 
                       table[, 6]))

           plot(brk,
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

           if(v)
             text(x=brk,
                  y=y,
                  labels=format(round(y,
                                      v.round), 
                                nsmall=v.round),
                  pos=v.pos, ...)
         })

  axis(1, 
       at=round(brk,
                x.round),
       las=xlas, ...)
}
