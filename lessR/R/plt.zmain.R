.plt.main <- 
function(x, y, by=NULL, n.cat=getOption("n.cat"),
         object="point", topic="data",

         col.fill=getOption("col.fill.pt"),
         col.stroke=getOption("col.stroke.pt"),
         col.bg=getOption("col.bg"),
         col.grid=getOption("col.grid"),
         col.box=getOption("col.box"),

         col.trans=NULL, col.area=NULL,

         cex.axis=0.76, col.axis="gray30", xy.ticks=TRUE,
         xlab=NULL, ylab=NULL, main=NULL, sub=NULL,
         value.labels=NULL, rotate.values=0, offset=0.5,

         size=NULL, shape="circle", means=TRUE, 
         sort.yx=FALSE, segments.y=FALSE, segments.x=FALSE,

         bubble.size=0.25, bubble.power=0.6, bubble.counts=TRUE,
         col.low=NULL, col.hi=NULL,

         fit.line="none", col.fit.line="gray55",

         ellipse=FALSE, col.ellipse="lightslategray",
         fill.ellipse="off", 

         method="overplot", pt.reg="circle", pt.out="circle", 
         col.out30="firebrick2", col.out15="firebrick4", new=TRUE,

         quiet=getOption("quiet"), fun.call=NULL, want.labels=TRUE, ...)  {


  # scale for regular R or RStudio
  adj <- .RSadj(bubble.size, cex.axis)
  bubble.size <- adj$bubble.size
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab
  size.txt <- adj$size.txt

  # want.labels set just for ttestPower, which provides its own labels
  # both x and y are plotted, even if only a single variable
  # for a 1-D bubble plot of a single factor var, y.call was set to 0's
  # numerical 1-D scatter plot done in .dp.main
  if (length(unique(y)) == 1)
    bubble1 <- TRUE
  else
    bubble1 <- FALSE

  unique.x <- ifelse(length(unique(x)) == length(x), TRUE, FALSE)
  unique.y <- ifelse(length(unique(y)) == length(y), TRUE, FALSE)

  # sort y by x option (intended for Cleveland dot plot)
  if (sort.yx) {
    if (!is.matrix(x))
      ord <- order(x)
    else
      if (ncol(x) == 2)
        ord <- order(x[,2] - x[,1])
      else {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Sorting not meaningful for more than two x-variables\n\n")
      }
    y <- factor(y, levels=y[ord])
  }

  do.ellipse <- ifelse(as.logical(ellipse[1]), TRUE, FALSE) 

  if (!is.null(value.labels)) value.labels <- gsub(" ", "\n", value.labels) 

  # all processing in terms of numeric variables
  # convert factors to numeric, save levels
  # x is always a matrix
  x.lvl <- NULL; y.lvl <- NULL  # if null, then not factors
  if (is.factor(x)) {
    x.lvl <- levels(x)
    if (is.null(value.labels)) value.labels <- gsub(" ", "\n", x.lvl) 
    x <- as.matrix(as.integer(x))
  }
  else {
    was.mat <- ifelse(is.matrix(x), TRUE, FALSE)
    x <- as.matrix(x)
    if (!was.mat) colnames(x) <- getOption("xname")
  }

  if (is.factor(y)) {
    y.lvl <- levels(y)
    y <- as.matrix(as.integer(y))
  }
  else {
    was.mat <- ifelse(is.matrix(y), TRUE, FALSE)
    y <- as.matrix(y)
    if (!was.mat) colnames(y) <- getOption("yname")
    #if (getOption("xname") == "Index") ylab <- ""
  }


  # see if trans is customized for this analysis
  if (is.null(col.trans)) {  # no change, so no trans for Cleveland dp
    if (is.null(x.lvl) && !is.null(y.lvl) && unique.y) { 
      trans.pts <- 0 
      col.fill <- .maketrans(col.fill, (1-trans.pts)*256)
    }
  }
  else {  # trans has been changed from default, so change col.fill
    trans.pts <- col.trans
    col.fill <- .maketrans(col.fill, (1-trans.pts)*256)
  }

  nrows <- length(x)
  # line chart
  size.pt <- NULL
  if (!is.null(xlab)) {
    if (xlab == "Index"  &&  object %in% c("line", "both")) {
      if (is.null(size)) size.pt <- .7
      if (is.null(size)) size.ln <- 1.25
    }
  }
  if (is.null(size.pt)) {
    sz <- ifelse (.Platform$OS == "windows", 1, 0.80)
    sz.pt <- ifelse (is.null(size), sz, size)
    size.pt <- ifelse (is.null(size), sz.pt, size)
    if (options("device") == "RStudioGD")
      size.pt <- ifelse (.Platform$OS == "windows", size.pt*1.45, size.pt*1.13)
    sz.ln <- ifelse (is.null(size), 2, size)
    size.ln <- ifelse (is.null(size), sz.ln, size)
  }

  if (is.null(col.fill)) col.fill <- "transparent"
  if (is.null(col.area)) col.area <- "transparent"

  num.cat.x <- is.null(x.lvl) && .is.num.cat(x, n.cat)
  cat.x <- ifelse (num.cat.x || !is.null(x.lvl), TRUE, FALSE)
  if (!bubble1) {
    num.cat.y <- is.null(y.lvl) && .is.num.cat(y, n.cat)
    cat.y <- ifelse (num.cat.y || !is.null(y.lvl), TRUE, FALSE)
  }
  else {
    num.cat.y <- FALSE
    cat.y <- FALSE
  }

  if (want.labels) {
    gl <- .getlabels(xlab, ylab, main, cex.lab=size.lab)
    x.name <- gl$xn; x.lbl <- gl$xl; x.lab <- gl$xb
    y.name <- gl$yn; y.lbl <- gl$yl; y.lab <- gl$yb
    main.lab <- gl$mb
    sub.lab <- gl$sb
    size.lab <- gl$cex.lab
    by.name <- getOption("byname")
  }
  else {
    x.lab <- xlab
    y.lab <- ylab
    main.lab <- main
    sub.lab <- sub
    size.lab <- 0.85
  }

  if (!is.null(x.name)) if (x.name == "Index") {
    if (ncol(y) > 1) y.lab <- ""
    if (!is.null(x.lbl)) y.lab <- paste(x.name, ": ", x.lbl, sep="")
  }

  # decimal digits
  digits.d <- .max.dd(y) + 1
  options(digits.d=digits.d)


  # -------------------------
  # plot
  # -------------------------
  # -------------------------

  x.val <- NULL
  y.val <- y.lvl  # if not reset to x value labels
  max.lbl.y <- NULL
  if (!is.null(value.labels)) { 
    x.val <- value.labels
    if (length(unique(y)) > 1) {  # see if set y axis values to those of x
      if (length(unique(na.omit(x))) == length(unique(na.omit(y)))) {
        if (all(sort(unique(x)) == sort(unique(y)))) {
          y.val <- value.labels
          v <- unlist(strsplit(value.labels, "\n", fixed=TRUE))
          max.lbl.y <- max(nchar(v))
        }
      }
    }
  }
  else {
    x.val <- x.lvl  # x.val is NULL if x is numeric, ignored
    y.val <- y.lvl  # y.val ignored if y is numeric 
  }


  # -----------------------
  # activate graphic system
  # -----------------------

  plot.new()  # activate for strwidth

  if (!is.null(y.val)) {
    max.width <- 0
    for (i in (1:length(y.val))) {
      li <- ifelse(!is.na(y.val[i]), strwidth(y.val[i], units="inches"), 0)
      if (li > max.width)
        max.width <- strwidth(y.val[i], units="inches")
    }
  }
  else 
    max.width <- strwidth(as.character(max(pretty(y))), units="inches")

  lm <- max(max.width + 0.3, 0.67)
  if (!is.null(y.lab)) {
    if (!nzchar(y.lab)[1]) lm <- lm - 0.3
    if (grepl("\n", y.lab[1], fixed=TRUE)) lm <- lm + .15
  }
  else
    lm <- lm - 0.3

  rm <- 0.25 
  if (!is.null(by)) rm <- 0.95
  if (ncol(y) > 1) rm <- 0.95
  tm <- ifelse (is.null(main), 0.25, 0.75)
  bm <- 0.75

  orig.params <- par(no.readonly=TRUE)
  on.exit(par(orig.params))
  RSpad.bm <- ifelse (options("device") == "RStudioGD", 0.4, 0)
  RSpad.m <- ifelse (options("device") == "RStudioGD", 0.3, 0)
  par(mai=c(bm+RSpad.bm, lm+RSpad.m, tm+RSpad.m, rm+RSpad.m))


  # -----------------------
  # setup coordinate system only with plot and type="n"
  # non-graphical parameters in ... generate warnings when no plot
  # -----------------------

  mn.x <- ifelse(is.null(x.lvl), min(x, na.rm=TRUE), 1)
  mx.x <- ifelse(is.null(x.lvl), max(x, na.rm=TRUE), length(x.lvl))
  mn.y <- ifelse(is.null(y.lvl), min(y, na.rm=TRUE), 1)
  mx.y <- ifelse(is.null(y.lvl), max(y, na.rm=TRUE), length(y.lvl))
  
  if (!do.ellipse) {

    if (cat.x) {
      mn.x <- mn.x - .5
      mx.x <- mx.x + .5
    }
    if (cat.y) {
      mn.y <- mn.y - .5
      mx.y <- mx.y + .5
    }

    if (topic %in% c("count", "prop")) mn.y <- 0
    if (topic == "prop") mx.y <- 1.0
    if (topic != "data" && (!all(y == 0))) mx.y <- mx.y + (.08 * (mx.y-mn.y))

    region <- matrix(c(mn.x, mx.x, mn.y, mx.y), nrow=2, ncol=2)

  }  # no ellipse

  else {  # set plot with sufficient room for ellipse and data
    cxy <- cor(x[,1],y[,1], use="complete.obs")
    s.x <- sd(x[,1], na.rm=TRUE); s.y <- sd(y[,1], na.rm=TRUE)
    m.x <- mean(x[,1], na.rm=TRUE); m.y <- mean(y[,1], na.rm=TRUE)
    lvl <- max(ellipse)
    region <- ellipse(cxy, scale=c(s.x, s.y), centre=c(m.x, m.y), level=lvl)
    region <- rbind(region, c(mn.x, mn.y))
    region <- rbind(region, c(mx.x, mx.y))
  }

  plot(region, type="n", axes=FALSE, ann=FALSE, ...)
  rm(region)


  # ----------------------
  # set up plot background
  # ----------------------

  # axis ticks and values
  if (cat.x) {
    if (!is.null(x.lvl)) axT1 <- 1:length(x.lvl)   # mark category values
    if (num.cat.x) axT1 <- sort(unique(x))
    #if (num.cat.x) axT1 <- axTicks(1)
  }
  else
    axT1 <- axTicks(1)  # else numeric, so all the ticks

  if (cat.y) { 
    if (!is.null(y.lvl)) axT2 <- 1:length(y.lvl)
    if (num.cat.y) axT2 <- sort(unique(y))
    #if (num.cat.y) axT2 <- axTicks(2)
  }
  else
    axT2 <- axTicks(2) 

  if (xy.ticks) {
    if (!bubble1)
      .axes(x.val, y.val, axT1, axT2,
            par("usr")[1], par("usr")[3], size.axis, col.axis,
            rotate.values, offset=offset, ...)
    else  # 1-l scatter plot of categorical variable
      .axes(x.val, NULL, axT1, NULL,
            par("usr")[1], par("usr")[3], size.axis, col.axis, 
            rotate.values, offset=offset, ...)
  }
       
  # axis labels 
  if (is.null(max.lbl.y)) {  # could be set earlier when x.val = y.val
    if (!is.null(y.lvl)) {
      max.lbl.y <- max(nchar(y.lvl))
    }
    else
      max.lbl.y <- max(nchar(axTicks(2)))
  }

  if (bubble1) {
    y.lab <- ""
    max.lbl.y <- 0
  }

  if (ncol(x) > 1) x.lab <- NULL
  .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl.y, 
          xy.ticks, offset=offset, cex.lab=size.lab, ...) 

  # color plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border=col.box)

  # grid lines
  abline(v=axT1, col=col.grid, lwd=.5)
  if (!bubble1)
    abline(h=axT2, col=col.grid, lwd=.5)

  # fill area under curve
  col.border <- ifelse (object=="point", col.stroke[1], "transparent")
  if (!is.null(col.area)  &&  col.area != "transparent") 
    polygon(c(x[1],x,x[length(x)]), c(min(y),y,min(y)),
            col=col.area, border=col.border)


  # ---------------
  # plot the values
  # ---------------

  n.col <- 1
  if (ncol(x) > 1) n.col <- ncol(x) 
  if (ncol(y) > 1) n.col <- ncol(y) 

  # plot lines and/or points
  if (object %in% c("point", "line", "both", "off")) {
    trans.pts <- ifelse(is.null(col.trans), getOption("trans.fill.pt"), col.trans)

    n.by <- ifelse (is.null(by), 0, nlevels(by))

    n.patterns <- max(n.col, n.by)

    stroke <- character(length=0)
    fill <- character(length=0)
    if (n.patterns == 1) {
      stroke[1] <- col.stroke[1]
      fill[1] <- col.fill[1]
    }
    else if (n.patterns == 2) {
      stroke[1] <- col.stroke[1]
      stroke[2] <- ifelse (length(col.stroke) > 1, col.stroke[2], col.stroke[1])
      fill[1] <- col.fill[1]
      if (object %in% c("line","both")) {
        fill[2] <- .col.discrete(bright=TRUE)[2]
        stroke[2] <- fill[2]
      }
      else {
        if (length(col.fill) == 1)
          fill[2] <- "transparent"
        else
          fill[2] <- col.fill[2]
      }
    }  # n.patterns=2
    else  # n.patterns > 2
      stroke <- .col.discrete(bright=TRUE)[1:n.patterns]

    if (n.patterns > 2)
      for (i in 1:length(stroke)) fill[i] <- .maketrans(stroke[i], (1-trans.pts)*256)

    # lines
    if (object == "line" || object == "both") {
      if (ncol(x) == 1  &&  ncol(y) == 1)
          lines(as.numeric(x[,1]),y[,1], col=fill, lwd=size.ln, ...)
      if (ncol(y) > 1)
        for (i in 1:ncol(y))
          lines(as.numeric(x[,1]),y[,i], col=fill[i], lwd=size.ln, ...)
      if (ncol(x) > 1)
        for (i in 1:ncol(x))
          lines(as.numeric(x[,i]),y[,1], col=fill[i], lwd=size.ln, ...)

      if (ncol(x) > 1)  # horizontal legend, on x-axis
        .plt.legend(colnames(x), TRUE, stroke, fill, shape, usr)

      if (ncol(y) > 1)  # vertical legend, on y-axis
        .plt.legend(colnames(y), FALSE, stroke, fill, shape, usr)
    }


    # points
    if (object == "point" || object == "both") {

      if (is.null(by)) { 

        if (object != "off") {

          if (ncol(y) == 1)
            for (i in 1:ncol(x))
              points(x[,i],y, pch=shape, col=stroke[i], bg=fill[i],
                              cex=size.pt, ...)
          else
            for (i in 1:ncol(y))
              points(x[,1],y[,i], pch=shape, col=stroke[i], bg=fill[i],
                              cex=size.pt, ...)

          if (ncol(x) > 1)  # horizontal legend, on x-axis
            .plt.legend(colnames(x), TRUE, stroke, fill, shape, usr)
          if (ncol(y) > 1)  # vertical legend, on y-axis
            .plt.legend(colnames(y), FALSE, stroke, fill, shape, usr)


          if (segments.y) { 
            if (ncol(x) == 1) # line segments from points to axis
              segments(x0=min(pretty(x)), y0=y, x1=x, y1=y, lty=1, lwd=.75, col=col.stroke)
            else if (ncol(x) == 2)  # line segments between points
              segments(x0=x[,1], y0=y[,1], x1=x[,2], y1=y[,1], lty=1, lwd=.75, col=col.stroke)
          }

          if (!(topic %in% c("count", "prop"))) {
            if (segments.x)
              segments(y0=par("usr")[3], x0=x, y1=y, x1=x, lty=1, lwd=.75,
                       col=col.stroke)
          }
          else {
            if (segments.x) 
              if (ncol(x) == 1)
                 segments(y0=0, x0=x, y1=y, x1=x, lty=1, lwd=1, col=col.stroke)
          }

        }  # object == "off"

        if (means) {
          pch.avg <- ifelse(getOption("colors")!="gray", 21, 23)
          bck.g <- ifelse(getOption("colors")!="gray", "gray15", "gray30")
          if (grepl(".black", getOption("colors"), fixed=TRUE)) bck.g <- "gray85"

          m.lvl <- numeric(length = 0)

          # plot means for factor x, num y
          if (!is.null(x.lvl) && is.null(y.lvl) && !unique.x) {
            for (i in (1:length(x.lvl))) 
              m.lvl[i] <- mean(y[x==i], na.rm=TRUE)
            abline(h=m.lvl, col="gray50", lwd=.5)
            points(m.lvl, pch=pch.avg, bg=bck.g)
          }

          # plot means for num x, factor y
          if (is.null(x.lvl) && !is.null(y.lvl) && !unique.y) {
            for (i in (1:length(y.lvl))) 
              m.lvl[i] <- mean(x[y==i], na.rm=TRUE)
            abline(v=m.lvl, col="gray50", lwd=.5)
            points(m.lvl, 1:length(y.lvl), pch=pch.avg, bg=bck.g)
          }
        }  # means
      }  # null by

      else {  # by grouping variable

        clr <- character(length(n.by))
        clr.tr <- character(length(n.by))  # translucent version of clr
        if (length(stroke) == 1) 
          for (i in 1:n.by) clr[i] <- stroke  # all levels get same color
        else
          clr <- stroke

        shp <- integer(length(n.by))
        if (length(shape) == 1)
          for (i in 1:n.by) shp[i] <- shape
        else
          shp <- shape

        shape.dft <- c(21,23,22,24,25,7:14)  # shape defaults
        if (length(stroke)==1 && length(shape)==1)  #  color, shape 
          for (i in 1:n.by) shp[i] <- shape.dft[i]  # fill is default shapes

        for (i in 1:n.by) {
          x.lv <- subset(x, by==levels(by)[i])
          y.lv <- subset(y, by==levels(by)[i])
          points(x.lv, y.lv, pch=shp[i], col=clr[i], bg=fill[i], cex=size.pt, lwd=0.75, ...)
        }

        .plt.by.legend(levels(by), stroke, fill, shp, trans.pts, col.bg, usr)

        }  # end by

      }  # object is point or both
    }  # object is point, line, both


  else if (object %in% c("bubble", "sunflower")) {

    n.patterns <- 1  # for fit.line

    mytbl <- table(x, y)  # get the counts

    # colors
    if (is.null(col.low) ||  is.null(col.hi)) {
      clr <- col.fill
      clr.stroke <- col.stroke
    }
    else {
      color.palette <- colorRampPalette(c(col.low, col.hi))
      clr <- color.palette(nrow(mytbl))
      clr.stroke <- "gray70"
    }

    # melt the table to a data frame
    k <- 0
    xx <- integer(length=0)
    yy <- integer(length=0)
    count <- integer(length=0)
    for (i in 1:nrow(mytbl)) {
      for (j in 1:ncol(mytbl)) {
        if (mytbl[i,j] != 0) {  # 0 plots to a single pixel, so remove
          k <- k + 1
          count[k] <- mytbl[i,j]
          xx[k] <- as.numeric(rownames(mytbl)[i])  # rownames are factors
          yy[k] <- as.numeric(colnames(mytbl)[j])
        }
      }
    }
    cords <- data.frame(xx, yy, count)

    c <- cords$count
    if (object == "bubble")
      symbols(as.numeric(cords$xx), as.numeric(cords$yy),
            circles=c**bubble.power, inches=bubble.size,
            bg=clr, fg=clr.stroke, add=TRUE, ...)

    else if (object == "sunflower") {
      sunflowerplot(cords$xx, cords$yy, number=c, 
        seg.col=col.stroke, col=col.fill, cex.axis=size.axis, col.axis=col.axis,
        xlab=x.lab, ylab=y.lab, add=TRUE)
  }

    if (bubble.counts  &&  object == "bubble") { 
      max.c <- max(c, na.rm=TRUE)  # bubble is too small for count
      #min.bubble <- 0.25 * max.c  # radius, bubble.power=1
      #min.bubble <- 0.10 * max.c  # bubble.power=.6
      min.bubble <- (bubble.power/9) * max.c
      for (i in 1:length(c))
        if (!is.na(c[i])) if (c[i] <= min.bubble) c[i] <- NA
      text(cords$xx, cords$yy, c, cex=size.txt)
    }
  }  # end bubble/sunflower 


  # ellipse option
  if (do.ellipse) {
    for (i in 1:length(ellipse)) {
      e <- ellipse(cxy, scale=c(s.x, s.y), centre=c(m.x, m.y), level=ellipse[i])
      polygon(e, border=col.ellipse, col=fill.ellipse, lwd=1.5)
    }
  }


  # fit line option
  if (fit.line != "none") { 

    for (i in 1:n.patterns) {

      if (n.patterns == 1) {  # one plot
        x.lv <- x[,1]
        y.lv <- y[,1]
        clr <- col.fit.line
      }

      else {  # multiple

        if (!is.null(by)) {  # multiple by plots
          x.lv <- subset(x, by==levels(by)[i])
          y.lv <- subset(y, by==levels(by)[i])
        }

        if (n.col > 1) {  # multiple variable plots
          if (ncol(x) > 1) {
            x.lv <- x[,i]
            y.lv <- y[,1]
          }
          else {
            x.lv <- x[,1]
            y.lv <- y[,i]
          }
        }

        clr <- ifelse (length(stroke) == 1, stroke, stroke[i])
      }  # multiple

      ln.type <- "solid"
      if (n.patterns == 2 && i == 2)
        ln.type <- "dashed"

      ok <- is.finite(x.lv) & is.finite(y.lv)
      if (any(ok)) {
        x.ok <- x.lv[ok]
        y.ok <- y.lv[ok]
        ord <- order(x.ok)
        x.ord <- x.ok[ord]
        y.ord <- y.ok[ord] 

        if (fit.line == "loess") 
          lines(x.ord, fitted(loess(y.ord~x.ord, ...)), col=clr, lwd=1.5, lty=ln.type)

        if (fit.line == "ls") {
          if(!is.factor(x.lv)) {
            model <- lm(y.ord ~ x.ord)
            abline(model$coef, col=clr, lwd=1.5, lty=ln.type, xpd=FALSE)
          }
          else cat("\n>>> Note: Least squares line not permitted for a factor.\n")
        }
      }
    }  # ith pattern
  }  # fit.line



  # -----------
  # text output
  # -----------

  if (!quiet  &&  topic == "data") {

    # function call for suggestions
    fncl <- .fun.call.deparse(fun.call) 
    fncl <- gsub(")$", "", fncl)  # get function call less closing ) 
    fncl <- gsub(" = ", "=", fncl)

    # traditional two-var numeric var scatter plot
    if (!cat.x  &&  !cat.y  &&  object %in% c("point", "bubble")) {

      for (i in 1:n.col) {

        txsug <- ""
        if (getOption("suggest")  &&  i == 1) {
          fc <- ""
          if (!grepl("ellipse", fncl)  &&  n.col == 1)
            fc <- paste(fc, ", ellipse=TRUE", sep="")
          if (!grepl("fit.line", fncl)) fc <- paste(fc, ", fit.line=TRUE", sep="")
          if (grepl("fit.line=TRUE", fncl))
            fncl <- sub("fit.line=TRUE", "fit.line=\"ls\"", fncl)
          if (nzchar(fc)) {
            fc <- paste(fncl, fc, ") ", sep="")
            txsug <- paste(">>> Suggest: ", fc, "\n")
          }
          fc <- ""
          if (!grepl("size", fncl)) fc <- paste(fc, ", size=3", sep="")
          if (!grepl("color.bg", fncl)) fc <- paste(fc, ", color.bg=\"off\"", sep="")
          if (!grepl("color.grid", fncl)) fc <- paste(fc, ", color.grid=\"off\"", sep="")
          if (grepl("ellipse", fncl)) fncl <- .rm.arg.l("ellipse", fncl) 
          if (nzchar(fc)) {
            fc <- paste(fncl, fc, ") ", sep="")
            fc <- sub(",,", ",", fc, fixed=TRUE)  # hack
            txsug <- paste(txsug,"\n>>> Suggest: ", fc)
          }
        }

        if (ncol(x) > 1) {
          options(xname = colnames(x)[i])
          stuff <- .cr.main(x[,i], y[,1], brief=TRUE, ...) 
        }
        else {
          options(yname = colnames(y)[i])
          stuff <- .cr.main(x[,1], y[,i], brief=TRUE, ...) 
        }

        txbck <- stuff$txb
        txdsc <- stuff$txd
        txinf <- stuff$txi

        class(txsug) <- "out_piece"
        class(txbck) <- "out_piece"
        class(txdsc) <- "out_piece"
        class(txinf) <- "out_piece"

        if (nzchar(txsug))
          output <- list(out_suggest=txsug, out_background=txbck,
            out_describe=txdsc, out_inference=txinf,
            r=stuff$r, tvalue=stuff$tvalue, df=stuff$df, pvalue=stuff$pvalue,
            lb=stuff$lb, ub=stuff$ub)
        else
          output <- list(out_background=txbck,
            out_describe=txdsc, out_inference=txinf,
            r=stuff$r, tvalue=stuff$tvalue, df=stuff$df, pvalue=stuff$pvalue,
            lb=stuff$lb, ub=stuff$ub)


        class(output) <- "out_all"
        print(output)
      }
    }


    # categorical var with numeric var, means plot or bubble-1D plot
    else if ((cat.x && !cat.y && !unique.x) || (!cat.x && cat.y && !unique.y)) {

      txsug <- ""
      if (getOption("suggest")) {

        fc <- ""
        if (!grepl("means", fncl)) fc <- paste(fc, ", means=FALSE", sep="")
        if (nzchar(fc)) {
          fc <- paste(fncl, fc, ") ", sep="")
          txsug <- paste(">>> To not plot the means:", fc, "\n")
        }
        fc <- ""
        if (!grepl("topic", fncl)) {
          fc <- paste(fc, ", topic=\"mean\"", sep="")
          if (grepl("means", fncl)) fncl <- .rm.arg.l("means", fncl) 
        }
        if (nzchar(fc)) {
          fc <- paste(fncl, fc, ") ", sep="")
          fc <- sub(",,", ",", fc, fixed=TRUE)  # hack
          txsug <- paste(txsug,"\n>>> To only plot the means: ", fc)
        }
      }

      if (!bubble1) {
        if (cat.x && !cat.y) {
          if (!is.null(x.lvl))  # convert back to a factor if was one
            x.by <- factor(x, levels=1:length(x.lvl), labels=x.lvl)
          else
            x.by <- x
          options(yname = x.name)  # reverse order of x and y for .ss.numeric
          options(xname = y.name)
          stats <- .ss.numeric(y, by=x.by, digits.d=digits.d, brief=TRUE)
        }
        else if (!cat.x && cat.y) {
          if (!is.null(y.lvl))  # convert back to a factor if was one
            y.by <- factor(y, levels=1:length(y.lvl), labels=y.lvl)
          else
            y.by <- y
          stats <- .ss.numeric(x, by=y.by, digits.d=digits.d, brief=TRUE)
        }

        txout <- stats$tx

        class(txout) <- "out_piece"

        output <- list(out_suggest=txsug, out_txt=txout)
        class(output) <- "out_all"
        print(output)
      }

      else {  # 1-D bubble plot of a factor var, y just a constant

      txsug <- ""
      if (getOption("suggest")) {
        fc <- ""
        if (!grepl("color.low", fncl))
          fc <- paste(fc, ", color.low=\"lemonchiffon2\"", sep="")
        if (!grepl("color.hi", fncl))
          fc <- paste(fc, ", color.hi=\"lightsteelblue2\"", sep="")
        if (nzchar(fc)) {
          fc <- paste(fncl, fc, ") ", sep="")
          fc <- sub(",,", ",", fc, fixed=TRUE)  # hack
          txsug <- paste(">>> Suggest: ", fc)
        }
      }

        if (!is.null(x.lvl))
          x.by <- factor(x, levels=1:length(x.lvl), labels=x.lvl)
        else
          x.by <- factor(x)

        stats <- .ss.factor(x.by, by=NULL, brief=TRUE, digits.d=NULL,
                            x.name, y.name, x.lbl, y.lbl)

        txttl <- stats$title
        counts <- stats$count
        chi <- stats$chi

        class(txsug) <- "out_piece"
        class(txttl) <- "out_piece"
        class(counts) <- "out_piece"
        class(chi) <- "out_piece"
        output <- list(out_suggest=txsug, out_title=txttl, out_counts=counts, out_chi=chi)
        class(output) <- "out_all"
        print(output)      
      }
    }


    # Cleveland dot plot
    else if ((cat.x && !cat.y && unique.x) || (!cat.x && cat.y && unique.y)) { 

      txsug <- ""
      if (getOption("suggest")) {
        fc <- ""
        if (!grepl("sort.yx", fncl)) fc <- paste(fc, ", sort.yx=TRUE", sep="")
        if (!grepl("segments.y", fncl)) fc <- paste(fc, ", segments.y=TRUE", sep="")
        if (!grepl("color.bg", fncl)) fc <- paste(fc, ", color.bg=\"off\"", sep="")
        if (!grepl("color.grid", fncl)) fc <- paste(fc, ", color.grid=\"off\"", sep="")
        if (nzchar(fc)) {
          fncl <- .fun.call.deparse(fun.call) 
          fncl <- gsub(")$", "", fncl)  # get function call less closing )
          fncl <- gsub(" = ", "=", fncl)
          fc <- paste(fncl, fc, ") ", sep="")
          fc <- sub(",,", ",", fc, fixed=TRUE)  # hack
          txsug <- paste(">>> Suggest: ", fc)
        }
      }

      if (!is.null(y.lvl))
        # convert back to a factor if was one originally
        y.by <- factor(y, levels=1:length(y.lvl), labels=y.lvl)
      else
        y.by <- y

        txout <- ""
        for (i in 1:ncol(x)) {
          stats <- .ss.numeric(x[,i], digits.d=digits.d, brief=TRUE)
          txout[length(txout)+1] <- paste("---", colnames(x)[i], "---")
          for (j in 2:length(stats$tx)) txout[length(txout)+1] <- stats$tx[j]
          if (i < ncol(x)) {
            txout[length(txout)+1] <- ""
            txout[length(txout)+1] <- ""
          }
        }

        class(txsug) <- "out_piece"
        class(txout) <- "out_piece"
        if (nzchar(txsug))
          output <- list(out_suggest=txsug, out_txt=txout)
        else
          output <- list(out_txt=txout)
        class(output) <- "out_all"
        print(output)
    }


    # categorical x and y vars
    else if (cat.x  &&  cat.y) {

      txsug <- ""
      if (getOption("suggest")) {
        fc <- ""
        if (!grepl("color.trans", fncl))
          fc <- paste(fc, ", color.trans=.8", sep="")
        if (!grepl("color.stroke", fncl))
          fc <- paste(fc, ", color.stroke=\"gray40\"", sep="")
        if (!grepl("bubble.counts", fncl))
          fc <- paste(fc, ", bubble.counts=FALSE", sep="")
        if (!grepl("color.bg", fncl))
          fc <- paste(fc, ", color.bg=\"off\"", sep="")
        if (!grepl("color.grid", fncl))
          fc <- paste(fc, ", color.grid=\"off\"", sep="")
        if (nzchar(fc)) {
          fncl <- .fun.call.deparse(fun.call) 
          fncl <- gsub(")$", "", fncl)  # get function call less closing )
          fncl <- gsub(" = ", "=", fncl)
          fc <- paste(fncl, fc, ") ", sep="")
          fc <- sub(",,", ",", fc, fixed=TRUE)  # hack
          txsug <- paste(">>> Suggest: ", fc)
        }
      }

      if (!is.null(x.lvl))
        x.fac <- factor(x, levels=1:length(x.lvl), labels=x.lvl)
      else
        x.fac <- x[,1]
      if (!is.null(y.lvl))
        y.fac <- factor(y, levels=1:length(y.lvl), labels=y.lvl)
      else
        y.fac <- y

      stats <- .ss.factor(x.fac, y.fac, brief=TRUE, digits.d=NULL,
                          x.name, y.name, x.lbl, y.lbl)

      txttl <- stats$txttl
      txfrq <- stats$txfrq
      txXV <- stats$txXV

      class(txsug) <- "out_piece"
      class(txttl) <- "out_piece"
      class(txfrq) <- "out_piece"
      class(txXV) <- "out_piece"
      output <- list(out_suggest=txsug, out_title=txttl, out_text=txfrq, out_XV=txXV)
      class(output) <- "out_all"
      print(output)
    }

  cat("\n")
  }       


}  # end plt.main

