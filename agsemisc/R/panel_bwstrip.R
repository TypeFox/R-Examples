panel.bwstrip = function (x, y, groups, subscripts, pch, col, 
      box.ratio = 5, varwidth = FALSE, whiskerpos = 0.1, logbase,
      type="mean,mad,strip,N,grid",
      densityplot=expression(density(X,cut=2)), strip.limit=100, 
      seplines = NULL, N.label = "N=% \n",
      extend = TRUE,
      levels.fos = NULL, ...) {
    #--- set up constants:
    box.dot <- trellis.par.get("box.dot")
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    strip.symbol <- trellis.par.get("superpose.symbol")
    reference.line <- trellis.par.get("reference.line")
    #--- preprocess arguments:
    y = as.numeric(y)
    groups = if(missing(groups)) factor(rep(TRUE,length(x)))
             else groups[subscripts]
    groupN = length(levels(groups))
    decompose = function(arg, colormode) {
      # extract pch characters and plot symbols from comma-separated string:
      # decompose("1,w,13,2,02",F) -> list("1", "w", 13, "2", 2)
      # decompose("1,red,13,2,02",T) -> list(1, "red", 13, 2, 2)
      # other scalars and vectors are just turned into a list.
      sep = ","
      if (is.character(arg) && length(arg) == 1 &&
          length(grep(sep, arg)) > 0) {
        strs = strsplit(arg, sep)[[1]]
        pattern = ifelse(colormode, "^[0-9]+$", "^[0-9][0-9]+$")
        integer = grep.b(pattern, strs) # agsemisc misc.R
        result = as.list(rep(NA, length(strs)))
        result[integer] = as.numeric(strs[integer]) # plotsymbol/color numbers
        result[!integer] = strs[!integer] # plot chars/color names
      } else {
        result = as.list(arg)
      }
      result
    }
    if (missing(pch))
      pch = rep(strip.symbol$pch, length.out=groupN)
    pch = decompose(pch, FALSE)
    if (missing(col))
      col = rep(strip.symbol$col, length.out=groupN)
    col = decompose(col, TRUE)
    box.ratio = if (box.ratio==1) 5 else box.ratio  # 1 is supplied as default
    if (missing(logbase)) logbase = FALSE
    if (length(type) == 1) type = strsplit(type, ",")[[1]]
    meanplot = "mean" %in% type
    madplot = "mad" %in% type
    stripplot = "strip" %in% type
    strip.outliersonly = strip.limit == TRUE
    if (strip.outliersonly) strip.limit = Inf
    show.N = "N" %in% type
    grid = "grid" %in% type || "g" %in% type
    #--- do the plotting for all boxplots in this panel:
    maxn <- max(by(x, y, length))
    yscale <- grid::current.viewport()$yscale
    if (is.null(levels.fos)) 
      levels.fos <- floor(yscale[2]) - ceiling(yscale[1]) + 1
    lower <- ceiling(yscale[1])
    height <- box.ratio/(1 + box.ratio)
    xscale <- grid::current.viewport()$xscale
    #--- separator lines:
    if (!is.null(seplines)) {
      gp = do.call("gpar", reference.line)
      sapply(seplines, function(x) 
        grid.lines(x = unit(c(0,1), "npc"), 
          y = unit(c(x,x), "native"), gp = gp))
    }
    if (grid) panel.grid(h=0, v=-1, col="grey50", lty=3)
    if (levels.fos > 0) {
        for (i in 1:levels.fos) {
            yval <- i
            #--- compute values for plot:
            X = x[y == yval]
            N = length(X)
            if (logbase)
              X <- logbase^X  # de-logarithmize X
            if (whiskerpos <= 0.25) {
              q <- quantile(X, c(1, 1-whiskerpos, 0.75, 0.5, 0.25, 
                                 whiskerpos, 0))
            } else { # whiskers based on interquartile range
              q <- quantile(X, c(1, 1, 0.75, 0.5, 0.25, 
                                 0, 0))
              bps = boxplot.stats(X, whiskerpos, FALSE, FALSE)
              q[2] = bps$stats[5] # upper whisker
              q[6] = bps$stats[1] # lower whisker
            }
            m <- mean(X)
            stderr <- ifelse(N > 1, sqrt(var(X)/N), 0)
            stderr.low  <- m - stderr
            stderr.high <- m + stderr
            mad.value <- ifelse(N > 1, mad(X)/sqrt(N),
                                0) # median absolute deviation
            mad.low   <- q[4] - mad.value
            mad.high  <- q[4] + mad.value
            do.densityplot = ("density" %in% type) && N > 3
            if (logbase) {
              X <- log(X, base=logbase)   # re-logarithmize X
              q <- log(q, base=logbase)
              m <- log(m, base=logbase)
              stderr.low  <- log(stderr.low,  base=logbase)
              stderr.high <- log(stderr.high, base=logbase)
              mad.low  <- log(mad.low,  base=logbase)
              mad.high <- log(mad.high, base=logbase)
            }
            stats <- list(stats = q[2:6], n = N,
                          out = X[X < q[6] | X > q[2]])
            whisker.low  = q[6]
            whisker.high = q[2]
            if (stats$n > 0) {
              grid::pushViewport(grid::viewport(y = unit(yval, "native"), 
                height = unit((if (varwidth) sqrt(stats$n/maxn) else 1)
                              * height, "native"), 
                xscale = xscale))
              r.x <- (stats$stats[2] + stats$stats[4])/2
              r.w <- stats$stats[4] - stats$stats[2]
              #--- plot box:
              gp = gpar(lwd = if (do.densityplot) 1 
                              else box.rectangle$lwd,
                     lty = box.rectangle$lty,
                     fill = if (do.densityplot) "transparent"
                            else box.rectangle$fill,
                     col = box.rectangle$col)
              grid.rect(x = unit(r.x, "native"), 
                        width = unit(r.w, "native"), gp = gp)
              gp = gpar(col = box.dot$col, cex = box.dot$cex)
              grid.points(x = stats$stats[3], y = 0.5,
                          pch=box.dot$pch, gp = gp)
              #--- plot whiskers:
              if (whiskerpos != 0.25) {
                gp = gpar(col = box.umbrella$col, 
                          lwd = box.umbrella$lwd, lty = box.umbrella$lty)
                grid.lines(x = unit(stats$stats[1:2], "native"), 
                           y = unit(c(0.5, 0.5), "npc"), gp = gp)
                grid.lines(x = unit(stats$stats[4:5], "native"), 
                           y = unit(c(0.5, 0.5), "npc"), gp = gp)
                grid.lines(x = unit(rep(stats$stats[1], 2), 
                             "native"), y = unit(c(0.3, 0.7), "npc"), gp = gp)
                grid.lines(x = unit(rep(stats$stats[5], 2), 
                             "native"), y = unit(c(0.3, 0.7), "npc"), gp = gp)
              }
              #--- meanplot:
              if (meanplot) {
                gp = gpar(col="black", lty=3, lwd=1)
                grid.lines(
                  x = unit(c(stderr.low, stderr.high), "native"),
                  y = unit(c(0.75, 0.75), "npc"), gp = gp)
                gp = gpar(col="black", cex=0.7)
                grid.points(x = m, y = 0.75, pch="M", gp = gp)
              }
              #--- madplot:
              if (madplot) {
                gp = gpar(col="black", lty=1, lwd=1)
                grid.lines(
                  x = unit(c(mad.low, mad.high), "native"),
                  y = unit(c(0.5, 0.5), "npc"), gp = gp)
              }
              #--- stripplot:
              if (stripplot) {
                if (strip.limit >= N) {
                  do.plot = function(idx, pch, col) {
                    if (length(idx) == 0) return(NULL)
                    gp = gpar(col=col)
                    grid.points(x = x[idx], 
                                y = runif(length(idx), 0.2, 0.3), 
                                pch = pch,
                                size=unit(strip.symbol$cex[1], "char"), 
                                gp = gp)
                  }
                  for (grpI in 1:groupN) {
                    grpname = levels(groups)[grpI]
                    isoutlier = x < whisker.low | x > whisker.high
                    choose = y==yval & groups==grpname &
                             (isoutlier | !strip.outliersonly)
                    do.plot((1:length(x))[choose], pch[[grpI]], col[[grpI]])
                  }
                }
              }
              #--- densityplot:
              if (do.densityplot) {
                dp <- eval(densityplot)
                gp = gpar(col = box.rectangle$col,
                          lwd = box.rectangle$lwd, lty = 1)
                grid.lines(x = unit(dp$x, "native"),
                           y = unit(dp$y/max(dp$y), "npc"), gp = gp)
                gp = gpar(col = strip.symbol$col[1],
                          lwd = 1, lty = 2)
                grid.lines(x = unit(range(X), "native"),
                           y = unit(c(0, 0), "npc"),
                           gp = gp) # base line
              }
              #--- show.N:
              if (show.N) {
                label = gsub("%", length(x[y==yval]), N.label)
                grid.text(label=label,
                          x = unit(1, "npc"), y = unit(0, "npc"),
                          just="right")
              }
              #--- extend:
              if (is.function(extend))
                extend(X, yval, groups, subscripts)
              if (extend == TRUE) {
                sig <- 4  # number of significant digits
                cat(yval, ":", format(q[c(7,5,4,3,1)], digits=sig),
                    "mean", format(m, digits=sig))
                # show quartile ratio or interquartile range:
                if(q[7] >= 0 && q[5] > 0)
                  cat(" qr", format(q[3]/q[5], digits=sig))
                else cat(" iqr", format(q[3]-q[5], digits=sig))
                cat(" N", N, "\n")
              }
              #--- finish:
              grid::popViewport()
            }
        }
    }
}
