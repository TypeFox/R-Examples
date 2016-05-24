## ----------------------------------------------------------------------------
## plot method for objects of class 'performance'
## ----------------------------------------------------------------------------

.get.arglist <- function( fname, arglist ) {
     if (fname=='plot')
       return(.select.args(arglist,
                           union(names(formals(plot.default)), names(par()))))
     else if (fname=='plot.xy')
       return(.select.args(arglist,
                           union( names(formals(plot.xy)), names(par()))))
     else return( .select.prefix( arglist, fname) )
 }

.downsample <- function( perf, downsampling ) {
    for (i in 1:length(perf@alpha.values)) {
        if (downsampling < 1 && downsampling > 0)
          ind <- round(seq(1, length(perf@alpha.values[[i]]),
                           length=(length(perf@alpha.values[[i]]) *
                                   downsampling)))
        else if (downsampling > 1)
          ind <- round(seq(1, length(perf@alpha.values[[i]]),
                           length=downsampling))
        else ind <- 1:length(perf@alpha.values[[i]])
        perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind]
        perf@x.values[[i]] <- perf@x.values[[i]][ind]
        perf@y.values[[i]] <- perf@y.values[[i]][ind]
    }
    return(perf)
}

.plot.performance <-
  function(perf, ..., avg="none",
           spread.estimate="none", spread.scale=1, show.spread.at=c(),
           colorize=FALSE, colorize.palette=rev(rainbow(256,start=0, end=4/6)),
           colorkey=colorize, colorkey.relwidth=0.25, colorkey.pos="right",
           print.cutoffs.at=c(),
           cutoff.label.function=function(x) { round(x,2) },
           downsampling=0, add=FALSE) {

      arglist <- c(lapply( as.list(environment()), eval ), list(...) )

      if (length(perf@y.values) != length(perf@x.values)) {
          stop("Performance object cannot be plotted.")
      }
      if (is.null(perf@alpha.values) && (colorize==TRUE ||
                                         length(print.cutoffs.at)>0)) {
          stop(paste("Threshold coloring or labeling cannot be performed:",
                     "performance object has no threshold information."))
      }
      if ((avg=="vertical" || avg=="horizontal") &&
          (colorize==TRUE || length(print.cutoffs.at)>0)) {
          stop(paste("Threshold coloring or labeling is only well-defined for",
                     "'no' or 'threshold' averaging."))
      }
    
      if (downsampling >0 ) perf <- .downsample( perf, downsampling)

      ## for infinite cutoff, assign maximal finite cutoff + mean difference
      ## between adjacent cutoff pairs
      if (length(perf@alpha.values)!=0) perf@alpha.values <-
        lapply(perf@alpha.values,
               function(x) { isfin <- is.finite(x);
                             x[is.infinite(x)] <-
                               (max(x[isfin]) +
                                mean(abs(x[isfin][-1] -
                                         x[isfin][-length(x[isfin])])));
                             x } )
      ## remove samples with x or y not finite
      for (i in 1:length(perf@x.values)) {
          ind.bool <- (is.finite(perf@x.values[[i]]) &
                       is.finite(perf@y.values[[i]]))
          
          if (length(perf@alpha.values)>0)
            perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]
          
          perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
          perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
      }
      arglist <- .sarg( arglist, perf=perf)
    
      if (add==FALSE) do.call( ".performance.plot.canvas", arglist )

      if (avg=="none") do.call(".performance.plot.no.avg", arglist)  
      else if (avg=="vertical")
        do.call(".performance.plot.vertical.avg", arglist)
      else if (avg=="horizontal")
        do.call(".performance.plot.horizontal.avg", arglist)
      else if (avg=="threshold")
        do.call(".performance.plot.threshold.avg", arglist)
  }

## ---------------------------------------------------------------------------
## initializing plots and plotting a canvas
## (can be skipped using 'plot( ..., add=TRUE)'
## ---------------------------------------------------------------------------

.performance.plot.canvas <- function(perf, avg, ...) {

    arglist <- list(...)

    axis.names <- list(x=perf@x.name, y=perf@y.name)
    if (avg=="horizontal" || avg=="threshold")
      axis.names$x <- paste("Average", tolower(axis.names$x))
    if (avg=="vertical" || avg=="threshold")
      axis.names$y <- paste("Average", tolower(axis.names$y))
    arglist <- .farg(arglist, xlab=axis.names$x, ylab=axis.names$y)

    arglist <-
      .farg(arglist,
            xlim=c(min(unlist(perf@x.values)), max(unlist(perf@x.values))),
            ylim=c(min(unlist(perf@y.values)), max(unlist(perf@y.values))))
    
    do.call("plot", .sarg(.slice.run(.get.arglist('plot', arglist)),
                          x=0.5, y=0.5, type='n', axes=FALSE))
    do.call( "axis", .sarg(.slice.run(.get.arglist('xaxis', arglist)),
                           side=1))
    do.call( "axis", .sarg(.slice.run(.get.arglist('yaxis', arglist)),
                           side=2))

    if (.garg(arglist,'colorkey')==TRUE) {
        colors <- rev( .garg(arglist,'colorize.palette') )
        max.alpha <- max(unlist(perf@alpha.values))
        min.alpha <- min(unlist(perf@alpha.values))
        col.cutoffs <- rev(seq(min.alpha,max.alpha, length=length( colors )))

        if ( .garg(arglist,'colorkey.pos')=="right") {
            
            ## axis drawing (ticks + labels)
            ## The interval [min.alpha,max.alpha] needs to be mapped onto
            ## the interval [min.y,max.y], rather than onto the interval
            ## [ylim[1],ylim[2]] ! In the latter case, NAs could occur in
            ## approxfun below, because axTicks can be out of the ylim-range
            ## ('yxaxs': 4%region)
            max.y <- max(axTicks(4))
            min.y <- min(axTicks(4))
            alpha.ticks <- .garg( arglist, c("coloraxis.at"))
            if (length(alpha.ticks)==0)
              alpha.ticks <- approxfun(c(min.y, max.y),
                                       c(min.alpha, max.alpha)) ( axTicks(4))
            alpha2y <- approxfun(c(min(alpha.ticks), max(alpha.ticks)),
                                 c(min.y,max.y))
            arglist <-
              .sarg(arglist,
                    coloraxis.labels=.garg(arglist,
                      'cutoff.label.function')(alpha.ticks),
                    coloraxis.at=alpha2y(alpha.ticks))
            
            do.call("axis",
                    .sarg(.slice.run(.get.arglist('coloraxis', arglist)),
                          side=4))

            ## draw colorkey
            ## each entry in display.bool corresponds to one rectangle of
            ## the colorkey.
            ## Only rectangles within the alpha.ticks range are plotted.
            ## y.lower, y.upper, and colors, are the attributes of the visible
            ## rectangles (those for which display.bool=TRUE)
            display.bool <- (col.cutoffs >= min(alpha.ticks) &
                             col.cutoffs < max(alpha.ticks))
            y.lower <- alpha2y( col.cutoffs )[display.bool]
            colors <- colors[display.bool]
            if (length(y.lower>=2)) {
                y.width <- y.lower[2] - y.lower[1]
                y.upper <- y.lower + y.width
                x.left <- .garg(arglist,'xlim')[2] +
                  ((.garg(arglist,'xlim')[2] - .garg(arglist,'xlim')[1]) *
                   (1-.garg(arglist,'colorkey.relwidth'))*0.04)
                x.right <- .garg(arglist,'xlim')[2] +
                  (.garg(arglist,'xlim')[2] -.garg(arglist,'xlim')[1]) * 0.04
                rect(x.left, y.lower, x.right, y.upper,
                     col=colors, border=colors,xpd=NA)
            }
        } else if (.garg(arglist, 'colorkey.pos') == "top") {
            ## axis drawing (ticks + labels)
            max.x <- max(axTicks(3))
            min.x <- min(axTicks(3))
            alpha.ticks <- .garg( arglist, c("coloraxis.at"))
            if (length(alpha.ticks)==0) {
                alpha.ticks <- approxfun(c(min.x, max.x),
                                         c(min.alpha, max.alpha))(axTicks(3))
            }
            alpha2x <- approxfun(c( min(alpha.ticks), max(alpha.ticks)),
                                 c( min.x, max.x))
            arglist <- .sarg(arglist,
                             coloraxis.labels=.garg(arglist,
                               'cutoff.label.function')(alpha.ticks),
                             coloraxis.at= alpha2x(alpha.ticks)) 
            do.call("axis",
                    .sarg(.slice.run( .get.arglist('coloraxis', arglist)),
                          side=3))

            ## draw colorkey
            display.bool <- (col.cutoffs >= min(alpha.ticks) &
                             col.cutoffs < max(alpha.ticks))
            x.left <- alpha2x( col.cutoffs )[display.bool]
            colors <- colors[display.bool]
            if (length(x.left)>=2) {
                x.width <- x.left[2] - x.left[1]
                x.right <- x.left + x.width
                y.lower <- .garg(arglist,'ylim')[2] +
                  (.garg(arglist,'ylim')[2] - .garg(arglist,'ylim')[1]) *
                    (1-.garg(arglist,'colorkey.relwidth'))*0.04
                y.upper <- .garg(arglist,'ylim')[2] +
                  (.garg(arglist,'ylim')[2] - .garg(arglist,'ylim')[1]) * 0.04
                rect(x.left, y.lower, x.right, y.upper,
                     col=colors, border=colors, xpd=NA)
            }
        }
    }
    
    do.call( "box", .slice.run( .get.arglist( 'box', arglist)))
}

## ----------------------------------------------------------------------------
## plotting performance objects when no curve averaging is wanted
## ----------------------------------------------------------------------------

.performance.plot.no.avg <- function( perf, ... ) {

    arglist <- list(...)
    arglist <- .farg(arglist, type= 'l')
    
    if (.garg(arglist, 'colorize') == TRUE) {
        colors <- rev( .garg( arglist, 'colorize.palette') )
        max.alpha <- max(unlist(perf@alpha.values))
        min.alpha <- min(unlist(perf@alpha.values))
        col.cutoffs <- rev(seq(min.alpha,max.alpha, length=length(colors)+1))
        col.cutoffs <- col.cutoffs[2:length(col.cutoffs)]
    }
    
    for (i in 1:length(perf@x.values)) {
        if (.garg(arglist, 'colorize') == FALSE) {
            do.call("plot.xy",
                    .sarg(.slice.run(.get.arglist('plot.xy', arglist), i),
                          xy=(xy.coords(perf@x.values[[i]],
                                        perf@y.values[[i]]))))
        } else {
            for (j in 1:(length(perf@x.values[[i]])-1)) {
                segment.coloring <-
                  colors[min(which(col.cutoffs <= perf@alpha.values[[i]][j]))]
                do.call("plot.xy",
                        .sarg(.slice.run(.get.arglist('plot.xy', arglist), i),
                              xy=(xy.coords(perf@x.values[[i]][j:(j+1)],
                                            perf@y.values[[i]][j:(j+1)])),
                              col= segment.coloring))
            }
        }

        print.cutoffs.at <- .garg(arglist, 'print.cutoffs.at',i)
        if (! is.null(print.cutoffs.at)) {
            text.x <- approxfun(perf@alpha.values[[i]], perf@x.values[[i]],
                                rule=2, ties=mean)(print.cutoffs.at)
            text.y <- approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                                rule=2, ties=mean)(print.cutoffs.at)
            do.call("points",
                    .sarg(.slice.run(.get.arglist('points', arglist),i),
                          x= text.x,
                          y= text.y))
            do.call("text",
                    .farg(.slice.run( .get.arglist('text', arglist),i),
                          x= text.x,
                          y= text.y,
                          labels=(.garg(arglist,
                                        'cutoff.label.function',
                                        i)(print.cutoffs.at))))
        }
    }
}

## ----------------------------------------------------------------------------
## plotting performance objects when vertical curve averaging is wanted
## ----------------------------------------------------------------------------

.performance.plot.vertical.avg <- function( perf, ...) {
    arglist <- list(...)
    arglist <- .farg(arglist,
                     show.spread.at= (seq(min(unlist(perf@x.values)),
                                          max(unlist(perf@x.values)),
                                          length=11)))
    perf.avg <- perf
    x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                    length=max( sapply(perf@x.values, length)))
    for (i in 1:length(perf@y.values)) {
        perf.avg@y.values[[i]] <-
          approxfun(perf@x.values[[i]], perf@y.values[[i]],
                    ties=mean, rule=2)(x.values)
    }
    perf.avg@y.values <- list(rowMeans( data.frame( perf.avg@y.values )))
    perf.avg@x.values <- list(x.values)
    perf.avg@alpha.values <- list()

    ## y.values at show.spread.at (midpoint of error bars )
    show.spread.at.y.values <-
      lapply(as.list(1:length(perf@x.values)),
             function(i) {
                 approxfun(perf@x.values[[i]], perf@y.values[[i]],
                           rule=2,
                           ties=mean)( .garg(arglist, 'show.spread.at'))
             })

    show.spread.at.y.values <- as.matrix(data.frame(show.spread.at.y.values ))
    colnames(show.spread.at.y.values) <- c()
    ## now, show.spread.at.y.values[i,] contains the curve y values at the
    ## sampling x value .garg(arglist,'show.spread.at')[i]
    
    if (.garg(arglist, 'spread.estimate') == "stddev" ||
        .garg(arglist, 'spread.estimate') == "stderror") {
        bar.width <- apply(show.spread.at.y.values, 1, sd)
        if (.garg(arglist, 'spread.estimate') == "stderror") {
            bar.width <- bar.width / sqrt( ncol(show.spread.at.y.values) )
        }
        bar.width <- .garg(arglist, 'spread.scale') * bar.width

        suppressWarnings(do.call("plotCI",
                                 .farg(.sarg(.get.arglist( 'plotCI', arglist),
                                             x=.garg(arglist,
                                               'show.spread.at'),
                                             y=rowMeans(
                                               show.spread.at.y.values),
                                             uiw= bar.width,
                                             liw= bar.width,
                                             err= 'y',
                                             add= TRUE),
                                       gap= 0,
                                       type= 'n')))
    }
    
    if (.garg(arglist, 'spread.estimate') == "boxplot") {
        do.call("boxplot",
                .farg(.sarg(.get.arglist( 'boxplot', arglist),
                            x= data.frame(t(show.spread.at.y.values)),
                            at= .garg(arglist, 'show.spread.at'),
                            add= TRUE,
                            axes= FALSE),
                      boxwex= (1/(2*(length(.garg(arglist,
                                                  'show.spread.at')))))))
        do.call("points",
                .sarg(.get.arglist( 'points', arglist),
                      x= .garg(arglist, 'show.spread.at'),
                      y= rowMeans(show.spread.at.y.values)))
    }
    
    do.call( ".plot.performance", .sarg(arglist,
                                       perf= perf.avg,
                                       avg= 'none',
                                       add= TRUE))
}

## ----------------------------------------------------------------------------
## plotting performance objects when horizontal curve averaging is wanted
## ----------------------------------------------------------------------------

.performance.plot.horizontal.avg <- function( perf, ...) {
    arglist <- list(...)
    arglist <- .farg(arglist,
                     show.spread.at= seq(min(unlist(perf@y.values)),
                       max(unlist(perf@y.values)),
                       length=11))
    perf.avg <- perf
    y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                    length=max( sapply(perf@y.values, length)))
    for (i in 1:length(perf@x.values)) {
        perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                            perf@x.values[[i]],
                                            ties=mean, rule=2)(y.values)
    }
    perf.avg@x.values <- list(rowMeans( data.frame( perf.avg@x.values )))
    perf.avg@y.values <- list(y.values)
    perf.avg@alpha.values <- list()
    
    ## x.values at show.spread.at (midpoint of error bars )
    show.spread.at.x.values <-
      lapply(as.list(1:length(perf@y.values)),
             function(i) {
                 approxfun(perf@y.values[[i]],
                           perf@x.values[[i]],
                           rule=2, ties=mean)(.garg(arglist,'show.spread.at'))
             } )
    show.spread.at.x.values <- as.matrix(data.frame(show.spread.at.x.values))
    colnames(show.spread.at.x.values) <- c()
    ## now, show.spread.at.x.values[i,] contains the curve x values at the
    ## sampling y value .garg(arglist,'show.spread.at')[i]
    
    if (.garg(arglist,'spread.estimate') == 'stddev' ||
        .garg(arglist,'spread.estimate') == 'stderror') {
        bar.width <- apply(show.spread.at.x.values, 1, sd)
        if (.garg(arglist,'spread.estimate')== 'stderror') {
            bar.width <- bar.width / sqrt( ncol(show.spread.at.x.values) )
        }
        bar.width <- .garg(arglist,'spread.scale') * bar.width

        suppressWarnings(do.call("plotCI",
                                 .farg(.sarg(.get.arglist( 'plotCI', arglist),
                                             x= rowMeans(
                                               show.spread.at.x.values),
                                             y= .garg(arglist,
                                               'show.spread.at'),
                                             uiw= bar.width,
                                             liw= bar.width,
                                             err= 'x',
                                             add= TRUE),
                                       gap= 0,
                                       type= 'n')))
    }
    
    if (.garg(arglist,'spread.estimate') == "boxplot") {
        do.call("boxplot",
                .farg(.sarg(.get.arglist( 'boxplot', arglist),
                            x= data.frame(t(show.spread.at.x.values)),
                            at= .garg(arglist,'show.spread.at'),
                            add= TRUE,
                            axes= FALSE,
                            horizontal= TRUE),
                      boxwex= 1/(2*(length(.garg(arglist,'show.spread.at'))))))
        do.call("points", .sarg(.get.arglist( 'points', arglist),
                                x= rowMeans(show.spread.at.x.values),
                                y= .garg(arglist,'show.spread.at')))
    }
    
    do.call( ".plot.performance", .sarg(arglist,
                                        perf= perf.avg,
                                        avg= 'none',
                                        add= TRUE))
}

## ----------------------------------------------------------------------------
## plotting performance objects when threshold curve averaging is wanted
## ----------------------------------------------------------------------------

.performance.plot.threshold.avg <- function( perf, ...) {
    arglist <- list(...)
    arglist <- .farg(arglist,
                     show.spread.at= seq(min(unlist(perf@x.values)),
                       max(unlist(perf@x.values)),
                       length=11))

    perf.sampled <- perf
    alpha.values <- rev(seq(min(unlist(perf@alpha.values)),
                            max(unlist(perf@alpha.values)),
                            length=max( sapply(perf@alpha.values, length))))
    for (i in 1:length(perf.sampled@y.values)) {
        perf.sampled@x.values[[i]] <-
          approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                    rule=2, ties=mean)(alpha.values)
        perf.sampled@y.values[[i]] <-
          approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                    rule=2, ties=mean)(alpha.values)
    }

    ## compute average curve
    perf.avg <- perf.sampled
    perf.avg@x.values <- list( rowMeans( data.frame( perf.avg@x.values)))
    perf.avg@y.values <- list(rowMeans( data.frame( perf.avg@y.values)))
    perf.avg@alpha.values <- list( alpha.values )
    
    x.values.spread <-
      lapply(as.list(1:length(perf@x.values)),
             function(i) {
                 approxfun(perf@alpha.values[[i]], perf@x.values[[i]],
                           rule=2, ties=mean)(.garg(arglist,'show.spread.at'))
             } )
    x.values.spread <- as.matrix(data.frame( x.values.spread ))
    y.values.spread <-
      lapply(as.list(1:length(perf@y.values)),
             function(i) {
                 approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                           rule=2, ties=mean)(.garg(arglist,'show.spread.at'))
             } )
    y.values.spread <- as.matrix(data.frame( y.values.spread ))

    if (.garg(arglist,'spread.estimate')=="stddev" ||
        .garg(arglist,'spread.estimate')=="stderror") {

        x.bar.width <- apply(x.values.spread, 1, sd)
        y.bar.width <- apply(y.values.spread, 1, sd)
        if (.garg(arglist,'spread.estimate')=="stderror") {
            x.bar.width <- x.bar.width / sqrt( ncol(x.values.spread) )
            y.bar.width <- y.bar.width / sqrt( ncol(x.values.spread) )
        }
        x.bar.width <- .garg(arglist,'spread.scale') * x.bar.width
        y.bar.width <- .garg(arglist,'spread.scale') * y.bar.width

        suppressWarnings( do.call("plotCI",
                                  .farg(.sarg(.get.arglist( 'plotCI', arglist),
                                              x= rowMeans(x.values.spread),
                                              y= rowMeans(y.values.spread),
                                              uiw= x.bar.width,
                                              liw= x.bar.width,
                                              err= 'x',
                                              add= TRUE),
                                        gap= 0,
                                        type= 'n')))
        
        suppressWarnings( do.call("plotCI",
                                  .farg(.sarg(.get.arglist( 'plotCI', arglist),
                                              x= rowMeans(x.values.spread), 
                                              y= rowMeans(y.values.spread),
                                              uiw= y.bar.width, 
                                              liw= y.bar.width, 
                                              err= 'y', 
                                              add= TRUE),
                                        gap= 0,
                                        type= 'n')))
    }

    if (.garg(arglist,'spread.estimate')=="boxplot") {
        do.call("boxplot",
                .farg(.sarg(.get.arglist('boxplot', arglist),
                            x= data.frame(t(x.values.spread)),
                            at= rowMeans(y.values.spread),
                            add= TRUE,
                            axes= FALSE,
                            horizontal= TRUE),
                      boxwex= 1/(2*(length(.garg(arglist,'show.spread.at'))))))
        do.call("boxplot",
                .farg(.sarg(.get.arglist('boxplot', arglist),
                            x= data.frame(t(y.values.spread)),
                            at= rowMeans(x.values.spread),
                            add= TRUE,
                            axes= FALSE),
                      boxwex= 1/(2*(length(.garg(arglist,'show.spread.at'))))))
        do.call("points", .sarg(.get.arglist('points', arglist),
                                x= rowMeans(x.values.spread),
                                y= rowMeans(y.values.spread)))
    }
    
    do.call( ".plot.performance", .sarg(arglist,
                                       perf= perf.avg,
                                       avg= 'none',
                                       add= TRUE))
}

