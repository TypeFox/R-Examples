##===================================================================
## Author : Y. Deville
##
## Return level plot from original data
##
## Should be transformed in a method valid for Renouv and Rendata
## objects by uncommenting the first part of the file, and changing
## the formal arguments of the second part.
##===================================================================

## RLplot <- function(x, ...) {
##   UseMethod("RLplot")
## }

## RLplot.Renouv <- function(x,
##                           pct.conf = x$pct.conf,
##                           show = list(OT = TRUE, quant = TRUE, conf = TRUE, MAX = TRUE,
##                               OTS = TRUE),
##                           mono = TRUE,
##                           predict = FALSE, 
##                           par = NULL,
##                           legend = TRUE, label = NULL,
##                           problim = NULL, Tlim = NULL,
##                           main = NULL, xlab = "periods", ylab = "level",
##                           posOptions = NULL,
##                           maxBlocks = 10L,
##                           ...) {

##     if (is.null(label)) label <- deparse(substitute(x))
    
##     plot.Renouv(x,
##                 pct.conf = pct.conf,
##                 show = show,
##                 mono = mono,
##                 predict = predict, 
##                 par = par,
##                 legend = legend, label = label,
##                 problim = problim, Tlim = Tlim,
##                 main = main, xlab = xlab, ylab = ylab,
##                 posOptions = posOptions,
##                 maxBlocks = maxBlocks,
##                 ...) 
       
## }

## RLplot.Rendata <- function(x,
##                            showHist = TRUE,
##                            mono = TRUE,
##                            par = NULL,
##                            legend = TRUE, label = NULL,
##                            Tlim = NULL,
##                            main = NULL, xlab = "periods", ylab = "level",
##                            posOptions = NULL,
##                            maxBlocks = 10L,
##                            ...) {
     
##     mc <- match.call(expand.dots = TRUE)
##     if (is.null(label)) label <- deparse(substitute(x))

##     if (maxBlocks > 10L) {
##         warning("'maxBlocks' must be <= 10")
##         maxBlocks <- 10L
##     }
    
##     if (!is.null(x$OTdata)) x$x.OT <- x$OTdata[ , x$info$varName]
##     else x$x.OT <- NULL
    
##     x$history.MAX <- makeMAXdata(x)
##     x$history.OTS <- makeOTSdata(x)
    
##     if (showHist) {
##         show = list(quant = FALSE, conf = FALSE, OT = TRUE, MAX = TRUE,
##             OTS = TRUE)
##     } else {
##         show = list(quant = FALSE, conf = FALSE, OT = TRUE, MAX = FALSE,
##             OTS = FALSE)
##     }
    
##     yLim <- range(x$x.OT, na.rm = TRUE)
    
##     if (showHist) {
##         yLim <- range(yLim, unlist(x$history.MAX$data),
##                       unlist(x$history.OTS$data))
##     }
##     labs <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
                  
##     if (is.null(main)) main <- ""

##     if (!is.null(Tlim)) { 
##         if ( !is.numeric(Tlim) || length(Tlim) != 2 || any(is.na(Tlim)) ||
##             any(Tlim < 0) ) {
##             stop("invalid limits in 'Tlim'.")
##         }
##         if (Tlim[1] < 0.01) Tlim[1] <- 0.01
##         xLim <-  log(Tlim)
##     } else xLim <- log(c(0.10, 500))
    
##     ## prepare (empty) plot
##     plot(x = xLim,
##          y = yLim,
##          type = "n",
##          main = main,
##          xlab = xlab,
##          ylab = ylab,
##          xaxt = "n",
##          ...)
##     ## add x-axis and grid
##     axis(side = 1, at = log(labs), labels = labs)
##     abline(v = log(labs), col = "gray", lty = "dotted")
##     abline(h = pretty(par()$usr[3:4]), col = "gray", lty = "dotted")
##     ## prepare legend

##     if (legend) RLlegend.ini()
  
##     ## Call 'lines' method
##     lines.Renouv(x = x,
##                  pct.conf = NULL,
##                  show = show,
##                  mono = mono,
##                  predict = FALSE,
##                  par = par,
##                  legend = legend,
##                  label = label,
##                  posOptions = posOptions,
##                  maxBlocks = maxBlocks,
##                  ...)
##     ## draw legend
##     if (legend) {
##         RLlegend.show()
##     }
## }


RLplot  <-  function(data,
                     x = NULL,
                     duration = 1,
                     lambda,
                     conf.pct = 95,
                     mono = TRUE,
                     mark.rl = 100,
                     mark.labels = mark.rl,
                     mark.col = NULL,
                     main = NULL,
                     ylim = NULL,  
                     ...) {
    if (mono) {
        l.cols <- c("black", "black")
        p.cols <- "black"
        l.typs <- c("solid", "dashed", "dotted")
    } else {
        l.cols <- c("SteelBlue4", "orangered", "SpringGreen3", "purple", "firebrick")
        p.cols <- "black"
        l.typs <- c("solid", "solid", "solid")
    }
    
    l.cols <- rep(l.cols, length.out = length(conf.pct)+1)
    l.typs <- rep(l.typs, length.out = length(conf.pct)+1)
    
    cnames <- colnames(data)
    candLnames <- match(paste("L", conf.pct, sep = "."), cnames)
    candUnames <- match(paste("U", conf.pct, sep = "."), cnames)
    
    nx <- length(x)
    
    ## freq.g = inverse return period (= freqency...)
    
    freq.g <- lambda*(1 - data$prob)
    x.g <- -log(freq.g)
    
    labs <- c(1, 2,  5, 10,20,  50, 100, 200, 500, 1000, 2000, 5000, 10000)
    
    ry <- c(min(data$L.95, na.rm = TRUE), max(data$U.95, na.rm = FALSE))
    
    if (!is.null(x)) {
        ry2 <- range(x, na.rm = TRUE)
        if (ry2[1] < ry[1]) ry[1] <- ry2[1]
        if (ry2[2] > ry[2]) ry[2] <- ry2[2]
    } 
    
    if (is.null(ylim)) ylim <- ry
    if (is.null(main)) main <- ""
    
    ## prepare plot
    plot(x = x.g,
         y = data$quant,
         type = "l",
         lwd = 2,
         main = main, ylim = ylim,
         xlab = "periods", ylab = "level",
         xaxt = "n",
         pch = 21,
         col = l.cols[1],
         ...)
    
    ind <- !is.na(candLnames) & !is.na(candUnames)
    
    for (i in 1:length(conf.pct))  {
        
        if (ind[i]) {
            iL <- candLnames[i]
            iU <- candUnames[i]
            lines(x = x.g, y = data[, iL],
                  type = "l", lwd = 2, col = l.cols[1+i], lty = l.typs[1+i])      
            lines(x = x.g, y = data[ ,iU],
                  type = "l", lwd = 2, col = l.cols[1+i], lty = l.typs[1+i])
        } else {
            warning("confidence limits for level ", conf.pct[i], "% not found in data")   
        }
    }   
        
    f.emp <- (1 - (1:nx)/(nx + 1)) * nx / duration
    
    ## abline(h = threshold, col = "SeaGreen3", lwd = 2)
    ## add empirical points if any
    
    if (!is.null(x)) {
        
        xs <- sort(x)
        
        points(x = -log(f.emp),
               y = xs,
               pch = 16, col = p.cols[1])
        
    }
    
    ## add x-axis
    axis(side = 1, at = log(labs), labels = labs)
    
    ## add upper x axis ??
    ## axis(side = 3, at = x.g, labels = data$prob)
    
    abline(v = log(labs), col = "gray", lty = "dotted")
    abline(h = pretty(ry), col = "gray", lty = "dotted")
    
    ## Not very useful...
    ##
    ## if (length(spec.x)) {
    ##   spec.f <- log(spec.rl)
    ##  points(x = spec.f,
    ##         y = spec.x,
    ##         pch = spec.pch,
    ##         bg = spec.bg,
    ##         cex = spec.cex,
    ##         col = spec.col)
    ## }
    
    Cols <- c("SteelBlue3", "orangered")
    
    if (sum(ind)) {
        legend(x = x.g[1],
               y = ry[2] -(ry[2] - ry[2])/10,
               legend = c("theo.", paste(conf.pct[ind], "%", sep = "")),
               col = l.cols[1:(1+sum(ind))],
               lwd = 2,
               lty = l.typs[1:(1+sum(ind))])
    } else {
        legend(x = x.g[1],
               y = ry[2] -(ry[2] - ry[2])/10,
               legend = "theo.",
               col = l.cols[1],
               lwd = 2,
               lty = l.typs[1])
    }
    
    if (!is.null(mark.rl)) {
        if ( mono || (is.null(mark.col)) ) mark.col <- rep("black", length(mark.rl))
        else  mark.col <- rep(mark.col, mength.out = length(mark.rl))
        
        for (i in 1:length(mark.rl)) {
            abline(v = log(mark.rl[i]), col = mark.col[i])
            text(x = log(mark.rl[i]),
                 y = ry[1] + (ry[2]-ry[1])/6,
                 labels = mark.labels[i],
                 col = mark.col[i],
                 pos = 4)
        }
    }
    
}

