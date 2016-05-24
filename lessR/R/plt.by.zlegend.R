# plot legend for by variables

.plt.by.legend <-
function(mylevels, stroke, fill, shp, trans.pts, col.bg, usr) {

  par(xpd=NA)  # allow drawing outside of plot region

  n.levels <- length(mylevels)
  by.name <- getOption("byname")

  legend.labels <- abbreviate(mylevels,6)
  legend.title  <- abbreviate(by.name, 12)
  ll <- legend(0,0, legend=legend.labels, title=legend.title, cex=.7,
           fill="gray", plot=FALSE)
  size <- (par("cxy")/par("cin"))  # 1 inch in user coordinates 
  epsilon <- (size[1] - ll$rect$w) / 2

  axis.vert <- usr[4] - usr[3]
  xleft <- usr[2] + epsilon   # usr[2] user coordinate of right axis
  lgnd.vhalf <- (ll$rect$h) / 2
  axis.cntr <- axis.vert / 2  + usr[3]
  ytop <- axis.cntr + lgnd.vhalf  # user coordinate of legend top

  if (trans.pts > 0.85) {  # points too light, reduce legend transparency
    legend.fill <- integer(length=n.levels)
    for (i in 1:n.levels) legend.fill[i] <- .maketrans(stroke[i],.7)
  }
  else 
    legend.fill <- fill

  colors <- getOption("colors")
  the.clr <- ifelse(grepl(".black", colors), "gray90", "black")

  legend(xleft, ytop, legend=legend.labels, title=legend.title, 
         pch=shp, horiz=FALSE, cex=.7, pt.cex=.9, box.lwd=.5, 
         box.col="gray30", bg=col.bg, col=stroke, pt.bg=fill, text.col=the.clr)

}
