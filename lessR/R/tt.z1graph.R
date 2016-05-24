.OneGraph <-
function(YA, bw1, Ynm, y.lbl, digits.d, brief, 
         n1, m1, mu0, mdiff, sw, smd, mmd, msmd,
         clpct, tvalue, pvalue, ub, lb,
         show.title, pdf.file, pdf.width, pdf.height) {

  dYA <- suppressWarnings(density(YA, bw1))

  if (!brief) {
    cat("\n\n------ Graphics Smoothing Parameter ------\n\n")
    mytitle <- "Density bandwidth for "
    cat(mytitle, .fmt(dYA$bw,digits.d), sep="", "\n")
  }


  cat("--------------------------------------------------\n")

  # values needed for graph
  min.x <- min(min(dYA$x), mu0)  # min x coordinate for graph
  max.x <- max(max(dYA$x), mu0)  # max x coordinate for graph
  max.y <- max(dYA$y)  # max y coordinate
  max.y <- max.y+.1*max.y  # allow room in graph region for d info

  # colors
  if (!grepl("gray", getOption("colors"))) {
    col.1 <- rgb(.63,.46,.15)
    col.m1 <- rgb(.71,.65,.65)
    col.1t <- rgb(.63,.46,.15, alpha=.7)
    col.1d <- rgb(.63,.46,.15, alpha=.4)
  }
  else {
    col.1 <- rgb(.40,.40,.40)
    col.m1 <- rgb(.40,.40,.40)
    col.1t <- rgb(.40,.40,.40, alpha=.7)
    col.1d <- rgb(.40,.40,.40, alpha=.4)
  }
  col.2 <- rgb(.49,.56,.69)

  # plot: set up coordinate system
  orig.params <- par(no.readonly=TRUE)
  on.exit(par(orig.params))
  par(mar=c(4.1,1.5,8,.4), mgp=c(3,.6,0), cex=.8, cex.axis=1, cex.lab=1.2)
  plot.new()
  plot.window(xlim=c(min.x,max.x), ylim=c(0,max.y))
  axis(1); box()
  if (nchar(y.lbl) > 50) y.lbl <- paste(substr(y.lbl,1,50), "...")
  title(xlab=y.lbl)

  xleft <- par("usr")[1]  # left side of graph
  xright <- par("usr")[2]  # right side of graph
  ybot <- par("usr")[3]  # bottom of graph
  ytop <- par("usr")[4]  # height of graph

  # vertical line for mean
  lines(c(m1,m1), c(0,ytop), lty="solid", lwd=.85, col=col.m1)
  if (!grepl("gray", getOption("colors")))
    lines(c(mu0,mu0), c(0,ytop), lty="twodash", lwd=.85, col=col.2)
  else
    lines(c(mu0,mu0), c(0,ytop), lty="twodash", lwd=.85, col=col.1)

  # curve area
  polygon(c(min(dYA$x),dYA$x,max(dYA$x)), c(0,dYA$y,0), col=col.1d, border=NA, 
          density=10, angle=45)

  # bottom border of density curve  
  segments(min(dYA$x), 0, max(dYA$x), 0, col=col.1)

  # density curve
  lwd.border <- 1.75
  if (.Platform$OS == "windows") lwd.border <- 2
  lines(dYA, col=col.1t, lty="solid", lwd=lwd.border)

  # minimum mean difference of practical importance
  if ( !is.null(mmd) | !is.null(msmd) ) {
    if (!is.null(mmd)) msmd <- mmd / sw
    if (!is.null(msmd)) mmd <- msmd * sw
    col.e <- "gray50"  # color for effect
    mid <- (m1 + mu0) / 2
    lr <- mid + .5*mmd  # line right
    ll <- mid - .5*mmd  # line left
    lines(c(lr,lr), c(ybot+.44*max.y,ytop-.44*max.y), lty="solid", lwd=2, col=col.e)
    lines(c(ll,ll), c(ybot+.44*max.y,ytop-.44*max.y), lty="solid", lwd=2, col=col.e)
    text(mid, ybot+.41*max.y, label=toString(.fmt(mmd,2)), col=col.e)
    text(mid, ytop-.41*max.y, label=toString(.fmt(msmd,2)), col=col.e)
    text(mid, ybot+.38*max.y, label="mmd", col=col.e)
    text(mid, ytop-.375*max.y, label="msmd", col=col.e)
  }

  # scale for s at top of graph
  mlow <- min(m1, mu0)
  mhi  <- max(m1, mu0)
  col.d.unit <- "gray35"
  # connect first seg to top
  segments(mlow, max.y-.01*max.y, mlow, ytop, lwd=1, col=col.d.unit) 
  # provide at least 2 labeled d units on sd scale at top
  max.i <- max(ceiling(abs(smd)), 2)
  for (i in 0:max.i) {  # sd scale at top
    x.i <- mlow+i*sw
    # sd units
    segments(x.i, max.y+.025*max.y, x.i, ytop, col=col.d.unit, lwd=1)
    # d units counted
    text(x.i, max.y+.01*max.y, labels=i)
    # horiz bar connects endpoints
    segments(mlow, ytop, x.i, ytop, col=col.d.unit, lwd=3)
    last.coord.x <- x.i
  }
  # connect last seg to top
  segments(last.coord.x, max.y+.025*max.y, last.coord.x, ytop, lwd=1, col=col.d.unit)
  # print d value towards top
  text((m1+mu0)/2, ytop-.07*max.y, label=.fmt(smd), col=col.d.unit, cex=.9)
  # horiz bar connects means
  segments(mlow, ytop-.09*max.y, mhi, ytop-.09*max.y, col=col.d.unit, lwd=1)
  # print d towards top
  text((m1+mu0)/2, ytop-.11*max.y, label="d", col=col.d.unit, cex=.9)

  # print mdiff value towards bottom  
  text((m1+mu0)/2, ybot+.11*max.y, label=.fmt(mdiff, digits.d), col=col.d.unit, cex=.9)
  # horiz bar connects means
  segments(mlow, ybot+.09*max.y, mhi, ybot+.09*max.y, col=col.d.unit, lwd=1)
  # print diff towards bottom
  text((m1+mu0)/2, ybot+.07*max.y, label="diff", col=col.d.unit, cex=.9)

  # title area, above graph
  if (show.title) {
    mtext(paste("One-Group Plot"), side=3, line=6.6, font=2)
    mtext(paste("Analyze",Ynm), side=3, line=5.4, font=3, cex=.8)
    mtext(bquote(paste("  t-test of mu0=", .(mu0), ":   t = ", .(.fmt(tvalue,3)), 
      ",  df = ", .(n1-1), ",   p-value = ", .(.fmt(pvalue,3)))), side=3, 
      line=3.9, cex=.8, adj=0)
    mtext(bquote(paste("  ",.(clpct), " Confidence Interval for Mean:  ",
      .(.fmt(lb,3)), " to ", .(.fmt(ub,3)))), side=3, line=2.8, cex=.8, adj=0)
    mtext(bquote(paste("  ", "n=", .(n1),
       "   m=", .(.fmt(m1, digits.d)),"   s=", .(.fmt(sw, digits.d)))),
       side=3, line=1.5, cex=.8, adj=0)
  }

}
