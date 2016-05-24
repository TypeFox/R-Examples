.dp.main <- 
function(x, by, size,
         col.fill, col.stroke, col.bg, col.grid, col.trans,
         shape.pts, cex.axis, col.axis, xlab, main, sub, cex,
         rotate.values, offset, method, pt.reg, pt.out, 
         col.out30, col.out15, quiet, new, ...) {


  # scale for regular R or RStudio
  adj <- .RSadj(cex.axis=cex.axis)
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab

  if (is.null(cex)) {
    sz <- ifelse (.Platform$OS == "windows", 1, 0.80)
    sz.pt <- ifelse (is.null(size), sz, size)
    size.pt <- ifelse (is.null(size), sz.pt, size)
    if (options("device") == "RStudioGD")
      size.pt <- ifelse (.Platform$OS == "windows", size.pt*1.45, size.pt*1.13)
  }

  if (is.factor(x)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The variable cannot be an R factor (categorical).\n")
  }

  # outlier points shapes and colors for gray scales
  if (getOption("colors") %in% c("gray", "gray.black")) {
    pt.out30 <- 23
    pt.out15 <- 22
  }
  else {
    pt.out30 <- pt.out
    pt.out15 <- pt.out
  }
  if (getOption("colors") == "gray") {
    col.out30 <- "black"
    col.out15 <- "gray30"
  }
  else if (getOption("colors") == "gray.black") {
    col.out30 <- "gray95"
    col.out15 <- "gray60"
  }

  # get values for ... parameter values
  stuff <- .getdots(...)
  col.main <- stuff$col.main
  col.lab <- stuff$col.lab
  col.sub <- stuff$col.sub
  cex.main <- stuff$cex.main

  # get variable labels if exist plus axes labels
  gl <- .getlabels(xlab, main=main, cex.lab=getOption("lab.size"))
  x.name <- gl$xn; x.lbl <- gl$xl
  x.lab <- gl$xb
  main.lab <- gl$mb
  sub.lab <- gl$sb
  cex.lab <- gl$cex.lab
  by.name <- getOption("byname")

  # text output (before remove missing)
  if (!quiet && new) .ss.numeric(x, brief=TRUE)

  n <- sum(!is.na(x))
  n.miss <- sum(is.na(x))
  if (n.miss > 0) x <- na.omit(x)
 
  if (new) {
    # set up plot area

    if (is.null(main) &&  is.null(by)) {
      orig.params <- par(no.readonly=TRUE)
      on.exit(par(orig.params))
      par(mar=c(4,4,2,2)+0.1)
    }

    if (!is.null(by)) {
      orig.params <- par(no.readonly=TRUE)
      on.exit(par(orig.params))
      par(omi=c(0,0,0,0.6))  # legend in right margin
    }

    stripchart(x, col="transparent", xlab=NULL, ylab=NULL, main=NULL,
       axes=FALSE, ann=FALSE, ...)

    # jitter passes to stripchart, but generates warning to axis
    #suppressWarnings(axis(1, cex.axis=cex.axis, col.axis=col.axis, ...))

    # axis, axis ticks
    .axes(x.lvl=NULL, y.lvl=NULL, axTicks(1), NULL,
          par("usr")[1], par("usr")[3], cex.axis, col.axis,
          rotate.values, offset, ...)

    # axis labels
    y.lab <- ""
    max.lbl <- max(nchar(axTicks(2)))
    .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl, 
          xy.ticks=TRUE, offset=offset, cex.lab=cex.lab, ...) 
    
    # colored background for plotting area
    usr <- par("usr")
    rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border="black")
    
    # grid lines computation and print
    vx <- pretty(c(usr[1],usr[2]))
    abline(v=seq(vx[1],vx[length(vx)],vx[2]-vx[1]), col=col.grid, lwd=.5)
  }

  # mark outliers
  q1 <- quantile(x, probs=0.25)
  q3 <- quantile(x, probs=0.75)
  lo30 <- q1 - 3.0*IQR(x)
  lo15 <- q1 - 1.5*IQR(x)
  up15 <- q3 + 1.5*IQR(x)
  up30 <- q3 + 3.0*IQR(x)
  stripchart(x[x<lo30], add=TRUE, method=method,
             col=col.out30, bg=col.out30, pch=pt.out30, cex=size, ...)
  stripchart(x[x<lo15], add=TRUE, method=method,
             col=col.out15, bg=col.out15, pch=pt.out15, cex=size, ...)
  stripchart(x[x>up15], add=TRUE, method=method,
             col=col.out15, bg=col.out15, pch=pt.out15, cex=size, ...)
  stripchart(x[x>up30], add=TRUE, method=method,
             col=col.out30, bg=col.out30, pch=pt.out30, cex=size, ...)

  # dp for regular points


  if (is.null(by)) {
    # see if trans is customized for this analysis
    if (!is.null(col.trans)) {
      trans.pts <- col.trans
      col.fill <- .maketrans(col.fill, (1-trans.pts)*256)
    }
    #trans.pts <- getOption("trans.fill.pt")
    #clr.trn <- .maketrans(col.fill, (1-trans.pts)*256)
    stripchart(x[x>lo15 & x<up15], add=TRUE, method=method,
                     col=col.stroke, pch=pt.reg, bg=col.fill, cex=size, ...)
  }

  else {  # by grouping variable
    n.levels <- nlevels(by)

    clr <- character(length(n.levels))
    if (length(col.stroke) == 1) 
      for (i in 1:n.levels) clr[i] <- col.stroke
    else
      clr <- col.stroke
    clr.tr <- clr

    shp <- integer(length(n.levels))
    if (length(shape.pts) == 1)
      for (i in 1:n.levels) shp[i] <- shape.pts
    else
       shp <- shape.pts
    shape.dft <- c(21,23,22,24,25,7:14)  # shape defaults
    if (length(col.stroke)==1 && length(shape.pts)==1)  # both shape and color default
      for (i in 1:n.levels) shp[i] <- shape.dft[i]  # fill with default shapes

    trans.pts <- getOption("trans.fill.pt")
    for (i in 1:n.levels) {
        clr.tr[i] <- .maketrans(clr.tr[i], (1-trans.pts)*256)
      x.lv <- subset(x, by==levels(by)[i])
      stripchart(x.lv, pch=shp[i], col=clr[i], bg=clr.tr[i], 
             cex=size.pt, lwd=0.75, add=TRUE, ...)
    }

    .plt.by.legend(levels(by), col.stroke, clr.tr, shp, trans.pts, col.bg, usr)

  }  # end by group

      txss <- ""
      if (!quiet) {
        digits.d <- NULL
        ssstuff <- .ss.numeric(x, digits.d=digits.d, brief=TRUE)
        txss <- ssstuff$tx

        txotl <- .outliers(x)
        if (length(txotl)==0) txotl <- "No outliers"

        class(txss) <- "out_piece"
        class(txotl) <- "out_piece"
        output <- list(out_ss=txss, out_outliers=txotl)
        class(output) <- "out_all"
        print(output)
      }

}
