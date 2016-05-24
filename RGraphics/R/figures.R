figure19.1 <- function() {
midpts <- barplot(1:10, col="gray90", axes=FALSE)
axis(2)
axis(1, at=midpts, labels=FALSE)


vps <- gridBase::baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.text(c("one", "two", "three", "four", "five",
            "six", "seven", "eight", "nine", "ten"), 
          x=unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=60)
popViewport(3)




}
figure19.2 <- function() {
hc <- hclust(dist(USArrests), "ave")
dend1 <- as.dendrogram(hc)
dend2 <- cut(dend1, h=70)



x <- 1:4
y <- 1:4
height <- factor(round(sapply(dend2$lower, 
                              attr, "height")))



space <- 1.2 * max(stringWidth(rownames(USArrests)))
dendpanel <- function(x, y, subscripts, ...) {
  pushViewport(viewport(gp=gpar(fontsize=8)),
               viewport(y=unit(0.95, "npc"), width=0.9,
                        height=unit(0.95, "npc") - space,
                        just="top"))
  par(plt=gridBase::gridPLT(), new=TRUE, ps=8)
  plot(dend2$lower[[subscripts]], axes=FALSE)
  popViewport(2)
}



trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
plot.new()
print(xyplot(y ~ x | height, subscripts=TRUE, 
             xlab="", ylab="",
             strip=strip.custom(style=4), 
             scales=list(draw=FALSE), 
             panel=dendpanel),
      newpage=FALSE)




}
figure3.1 <- function() {
par(oma=rep(3, 4), bg="gray80")
plot(c(0, 1), c(0, 1), type="n", ann=FALSE, axes=FALSE)
box("outer", col="gray")
# set clipping to figure region
par(xpd=TRUE)
# deliberately draw a stupidly large rectangle
rect(-1, -1, 2, 2, col="gray90")
box("figure")
# set clipping back to plot region
par(xpd=FALSE)
# deliberately draw a stupidly large rectangle
rect(-1, -1, 2, 2, col="gray80")
box("plot", lty="dashed")
text(.5, .5, "Plot Region")
mtext("Figure Region", side=3, line=2)
for (i in 1:4)
    mtext(paste("Outer margin", i), side=i, line=1, outer=TRUE)



}
figure3.2 <- function() {
par(oma=rep(3, 4), mfrow=c(3,2), bg="gray80")
for (i in 1:6) {
    if (i == 3) {
      omar <- par(mar=c(2, 2, 2, 1))  
      plot(c(0, 1), c(0, 1), type="n", ann=FALSE, axes=FALSE)
      par(xpd=TRUE)
      rect(-1, -1, 2, 2, col="gray90")
      box("figure")
      par(xpd=FALSE)
      rect(-1, -1, 2, 2, col="gray80")
      box("plot", lty="dashed")
      text(.5, .5, "Current Plot Region", cex=1.5)
      mtext("Current Figure Region", side=3)
      par(omar)
    } else {
      omar <- par(mar=rep(0, 4))  
      plot(c(0, 1), c(0, 1), type="n", ann=FALSE, axes=FALSE)
      par(xpd=TRUE)
      rect(-1, -1, 2, 2, col="gray90")
      box("figure")
      text(.5, .5, paste("Figure", i), cex=1.5)
      par(omar)
    }
}
box("outer", col="gray")
for (i in 1:4)
    mtext(paste("Outer margin", i), side=i, line=1, outer=TRUE)



}
figure3.3 <- function() {
par(mar=c(3, 6, 2, 2), xaxs="i", yaxs="i", xpd=FALSE, las=1)
    plot(c(0, 1), c(0, 1), type="n", ann=FALSE, axes=FALSE)
    box("figure")
    rect(0, 0, 1, 1, col="light gray", border="gray")
    axis(1, at=c(0, 1), c("", ""))
    mtext("Min x-value", side=1, adj=0, line=1)
    mtext("Max x-value", side=1, adj=1, line=1)
    axis(2, at=c(0, 1), c("", ""))
    mtext("Min y-value", side=2, at=0, adj=1, line=1)
    mtext("Max y-value", side=2, at=1, adj=1, line=1)
    lines(c(.6, .6, 0), c(0, .6, .6), lty="dashed")
    text(.6, .6, expression(paste("The location ", 
            group("(",list(x[i], y[i]),")"))), pos=3)
    points(.6, .6, pch=16)
    axis(1, at=.6, "")
    mtext(expression(x[i]), side=1, at=.6, line=.7)
    axis(2, at=.6, "")
    mtext(expression(y[i]), side=2, at=.6, line=.7)
        



}
figure3.4 <- function() {
pushViewport(viewport(layout=grid.layout(3, 1, 
  heights=unit(rep(1, 3), c("null", "cm", "null")))))
pushViewport(viewport(layout.pos.row=1))
grid.rect()
pushViewport(plotViewport(c(5, 5, 3, 2), xscale=c(0, 11)))
grid.rect(gp=gpar(col="gray"))
grid.text("Current Plot", gp=gpar(col="gray"))
grid.rect(0, unit(-5, "lines"), 1, unit(5, "lines"),
          just=c("left", "bottom"), gp=gpar(col="gray", fill="light gray"))
grid.text("Figure\nMargin\n1", y=unit(-2.5, "lines"))
grid.lines(c(0, 1), c(0, 0))
grid.segments(c(0, 1), c(0, 0), c(0, 1), unit(c(.5, .5), "lines"))
grid.text(c("xmin", "xmax"), c(0, 1), unit(c(1, 1), "lines"))
grid.lines(c(0, 0), unit(c(-1, -4), "lines"))
grid.segments(c(0, 0), unit(c(-1, -4), "lines"), 
              unit(c(-.5, -.5), "lines"), unit(c(-1, -4), "lines"))
grid.text(c("0 lines", "3 lines"),
          unit(c(-1, -1), "lines"), unit(c(-1, -4), "lines"),
          just=c("right", "bottom"))
popViewport(2)
pushViewport(viewport(layout.pos.row=3))
grid.rect()
pushViewport(plotViewport(c(5, 5, 3, 2), yscale=c(0, 11)))
grid.rect(gp=gpar(col="gray"))
grid.text("Current Plot", gp=gpar(col="gray"))
grid.rect(unit(-5, "lines"), 0, unit(5, "lines"), 1,
          just=c("left", "bottom"), gp=gpar(col="gray", fill="light gray"))
grid.text("Figure\nMargin\n2", x=unit(-2.5, "lines"))
grid.lines(c(0, 0), c(0, 1))
grid.segments(c(0, 0), c(0, 1), unit(c(.5, .5), "lines"), c(0, 1))
grid.text(c("ymin", "ymax"), unit(c(1, 1), "lines"), c(0, 1), just="left")
grid.lines(unit(c(0, -3), "lines"), c(0, 0))
grid.segments(unit(c(0, -3), "lines"), c(0, 0), 
              unit(c(0, -3), "lines"), unit(c(-.5, -.5), "lines"))
grid.text(c("0 lines", "3 lines"),
          unit(c(0, -3), "lines"), unit(c(-1, -1), "lines"),
          rot=90, just=c("right", "bottom"))
popViewport(2)
popViewport()



}
figure3.5 <- function() {
pushViewport(viewport(layout=grid.layout(3, 1, 
  heights=unit(c(1, 1, 1), c("null", "cm", "null")))))

# First page
pushViewport(viewport(layout.pos.row=3, 
  layout=grid.layout(3, 4, 
    widths=unit(c(2.5, 1, 1, 1), c("cm", "null", "null", "cm")),
    heights=unit(c(1, 1, 2.5), c("cm", "null", "cm")))))
grid.rect(gp=gpar(col="black"))
for (i in 2) {
  for (j in 2:3) {
    pushViewport(viewport(layout.pos.col=j, layout.pos.row=i))
    grid.rect(gp=gpar(col="gray"))
      pushViewport(plotViewport(c(2, 2, 1, 1), xscale=c(0, 11),
        gp=gpar(col="gray")))
      grid.rect(gp=gpar(col="gray"))
      grid.text(paste("Plot", j - 1))
      popViewport()      
    popViewport()
  }
}
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
grid.rect(gp=gpar(col="gray", fill="light gray"))
grid.text("Outer\nMargin\n2")
grid.lines(c(1, 1), c(0, 1))
grid.segments(c(1, 1), c(0, 1), 
              unit(1, "npc") + unit(c(.5, .5), "lines"),
              c(0, 1))
grid.text(0:1, 
          unit(1, "npc") + unit(c(1, 1), "lines"),
	  c(0, 1))
grid.lines(unit(1, "npc") - unit(c(0, 3), "lines"), c(0, 0))
grid.segments(unit(1, "npc") - unit(c(0, 3), "lines"), 
	      c(0, 0), 
              unit(1, "npc") - unit(c(0, 3), "lines"), 
	      unit(c(-.5, -.5), "lines"))
grid.text(c("0 lines", "3 lines"),
          unit(1, "npc") - unit(c(0, 3), "lines"),
          unit(c(-1, -1), "lines"),
          rot=90, just=c("right", "bottom"))
popViewport(2)

# Second page
pushViewport(viewport(layout.pos.row=1, 
  layout=grid.layout(3, 4, 
    widths=unit(c(2.5, 1, 1, 1), c("cm", "null", "null", "cm")),
    heights=unit(c(1, 1, 2.5), c("cm", "null", "cm")))))
grid.rect(gp=gpar(col="black"))
for (i in 2) {
  for (j in 2:3) {
    pushViewport(viewport(layout.pos.col=j, layout.pos.row=i))
    grid.rect(gp=gpar(col="gray"))
      pushViewport(plotViewport(c(2, 2, 1, 1), xscale=c(0, 11),
        gp=gpar(col="gray")))
      grid.rect(gp=gpar(col="gray"))
      grid.text(paste("Plot", j - 1))
      popViewport()      
    popViewport()
  }
}
pushViewport(viewport(layout.pos.row=3, layout.pos.col=2:3))
grid.rect(gp=gpar(col="gray", fill="light gray"))
grid.text("Outer Margin 1")
grid.lines(c(0, 1), c(1, 1))
grid.segments(c(0, 1), 
              unit(c(1, 1), "npc"),
	      c(0, 1),
              unit(c(1, 1), "npc") + unit(.5, "lines"))
# grid.rect(c(0, 1), 
# 	  unit(c(1, 1), "npc") + unit(1, "lines"),
# 	  unit(c(1, 1), "strwidth", list("0", "1")),
# 	  unit(c(1, 1), "strheight", list("0", "1")),
# 	  gp=gpar(col=NULL, fill="white"))
grid.text(c(0, 1), 
          c(0, 1),
          unit(c(1, 1), "npc") + unit(1, "lines"))
grid.lines(c(0, 0), unit(1, "npc") - unit(c(1, 4), "lines"))
grid.segments(c(0, 0), 
              unit(1, "npc") - unit(c(1, 4), "lines"), 
	      unit(c(-.5, -.5), "lines"),
              unit(1, "npc") - unit(c(1, 4), "lines"))
grid.text(c("0 lines", "3 lines"),
          unit(c(-1, -1), "lines"),
	  unit(1, "npc") - unit(c(1, 4), "lines"),
          just=c("right", "bottom"))
popViewport(2)

popViewport()



}
figure3.6 <- function() {
par(mar=rep(0, 4), cex=0.7)
plot.new()
plot.window(c(0.05, 0.95), 0:1)
family <- c("sans", "serif", "mono")
face <- 1:4
for (i in 1:4)
  for (j in 1:3) {
    par(family=family[j], lheight=1.5)
    text(seq(.15, .85, length=4)[i],
         seq(.25, .75, length=3)[j],
         paste("family=\"", family[j], "\"\nfont=", face[i], sep=""),
         font=face[i])
  }
segments(.02, c(.125, .375, .625, .875), 
         .98, c(.125, .375, .625, .875), col="gray")
segments(.02, c(.125, .375, .625, .875) - .01, 
         .02, c(.125, .375, .625, .875) + .01, col="gray")
segments(.98, c(.125, .375, .625, .875) - .01, 
         .98, c(.125, .375, .625, .875) + .01, col="gray")
rect(c(.27, .5, .73) - .01,
     .1,
     c(.27, .5, .73) + .01,
     .9, col="white", border=NA)



}
figure3.7 <- function() {
par(mar=rep(0, 4), xaxs="i", yaxs="i", cex=0.8)
plot.new()
par(new=TRUE)
grid.rect(gp=gpar(col="gray"))
ncol <- 4
nrow <- 4
xadj <- c(1, 0.5, NA, 0)
yadj <- c(1, 0.5, NA, 0)
size <- unit(3, "mm")
for (i in 1:nrow) {
  for (j in 1:ncol) {
    x <- i/(nrow + 1)
    y <- j/(ncol + 1)
    xu <- unit(x, "npc")
    yu <- unit(y, "npc")
    grid.segments(unit.c(xu - size, xu),
                  unit.c(yu, yu - size),
                  unit.c(xu + size, xu),
                  unit.c(yu, yu + size),
		  gp=gpar(col="gray"))
    text(x, y, paste("c(", xadj[j], ", ", yadj[i], ")", sep=""),
         adj=c(xadj[j], yadj[i]))
  }
}



}
figure3.8 <- function() {
ncol <- 6
nrow <- 1
grid.rect(gp=gpar(col="gray"))
for (i in 1:nrow) {
  for (j in 1:ncol) {
    x <- unit(j/(ncol+1), "npc")
    y <- unit(i/(nrow + 1), "npc")
    pch <- (i - 1)*ncol + j - 1
    grid.points(x + unit(3, "mm"), y, 
      pch=pch, gp=gpar(fill="gray"))
    grid.text(pch, x - unit(3, "mm"), y, gp=gpar(col="gray"))
  }
}



}
figure3.9 <- function() {
x <- -5:5
y <- -x^2 + 25
plottype <- function(type) {
  par(mar=c(1, 0, 1, 0), pty="s")
  plot.new()
  plot.window(c(-6, 6), c(-2, 27))
  box(col="gray")
  points(x, y, type=type)
  mtext(paste("type=\"", type, "\"", sep=""))
}



par(mfrow=c(3, 2))
plottype("p")
plottype("l")
plottype("b")
plottype("o")
plottype("h")
plottype("s")



}
figure3.10 <- function() {
axisfun <- function(mgp=c(3, 1, 0), xaxs="r", tcl=-.5,
                    mgpcol="black", xaxscol="black", tclcol="black") {
  par(mar=c(5, 1, 0, 1), mgp=mgp, xaxs=xaxs, tcl=tcl, pty="s")
  plot.new()
  box(col="gray")
  text(.5, .75, paste("mgp=c(", paste(mgp, collapse=", "), ")", sep=""),
       col=mgpcol)
  text(.5, .5, paste("xaxs=\"", xaxs, "\"", sep=""),
       col=xaxscol)
  text(.5, .25, paste("tcl=", tcl, sep=""),
       col=tclcol)
  axis(1, at=c(0, .5, 1))
  title(xlab="X-axis Label")
}



par(mfrow=c(2, 2))
axisfun()
axisfun(mgp=c(2, 0.3, 0), tcl=0.2, xaxscol="gray")
axisfun(xaxs="i", mgpcol="gray", tclcol="gray")



}
figure3.11 <- function() {

par(oma=rep(3, 4))
vps <- gridBase::baseViewports()
# Annotation helper function
annWidth <- function(x, y, lab, above=TRUE, horiz=TRUE) {
  grid.lines(x=x, y=y, 
             arrow=arrow(ends="both", angle=10, type="closed",
                         length=unit(3, "mm")), 
             gp=gpar(fill="black"))
  nl <- length(lab)
  if (nl > 1) {
    y <- y + unit(c(-0.5, 0.5), "lines")
    if (horiz) {
      vjust <- 1:0
      hjust <- 0.5
      rot <- 0
    } else {
      hjust <- 1:0
      vjust <- 0.5
      rot <- 90
    }
  } else {
    hjust <- 0.5
    rot <- 0
    if (above) {
      y <- y + unit(0.5, "lines")
      vjust <- 0
    } else {
      y <- y - unit(0.5, "lines")
      vjust <- 1
    }
  }
  grid.text(lab,
            x=0.5*sum(x),
            y=y, hjust=hjust, vjust=vjust, rot=rot,
            gp=gpar(fontfamily="mono", cex=1))
}
# Annotate whole page
grid.rect(gp=gpar(col="gray", fill="gray80"))
annWidth(0:1, unit(1, "npc") - unit(1.5, "lines"), "din[1]")
# grid.lines(x=0.5)
annWidth(unit(c(0, 3), "lines"), unit(0.7, "npc"), c("omi[2]", "oma[2]"))
annWidth(unit(1, "npc") - unit(c(0, 3), "lines"),
         unit(0.7, "npc"), c("omi[4]", "oma[4]"))
annWidth(unit(c(0, 3), "lines"), unit(0.3, "npc"), 
         "omd[1]", above=FALSE)
annWidth(unit.c(unit(0, "npc"),
                unit(1, "npc") - unit(3, "lines")),
         unit(2, "lines"), "omd[2]",
         above=FALSE)
# Annotate figure region
pushViewport(do.call("vpStack", vps[1:2]))
grid.rect(gp=gpar(fill="gray90"))
annWidth(0:1, unit(1, "npc") - unit(1.5, "lines"), "fin[1]")
annWidth(unit(c(0, 4.1), "lines"), unit(0.6, "npc"), c("mai[2]", "mar[2]"))
annWidth(unit(1, "npc") - unit(c(0, 2.1), "lines"),
         unit(0.6, "npc"), c("mai[4]", "mar[4]"), horiz=FALSE)
annWidth(unit(c(0, 4.1), "lines"), unit(0.4, "npc"), 
         "plt[1]", above=FALSE)
annWidth(unit.c(unit(0, "npc"),
                unit(1, "npc") - unit(2.1, "lines")),
         unit(4, "lines"), "plt[2]",
         above=FALSE)
# Annotate plot region
pushViewport(vps[[3]])
grid.rect(gp=gpar(lty="dashed", fill="gray80"))
annWidth(0:1, unit(1, "npc") - unit(1.5, "lines"), "pin[1]")
popViewport(3)



}
figure3.12 <- function() {
grid.lshow <- function(i, j, lab, order, nrow, ncol, heights, respect) {
  pushViewport(viewport(layout.pos.col=j, layout.pos.row=i))
  pushViewport(viewport(width=unit(1, "npc") - unit(2, "lines"),
               height=unit(1, "npc") - unit(3, "lines"),
	       y=unit(3, "lines"), just="bottom", 
    layout=grid.layout(nrow, ncol, heights=heights, 
      respect=respect)))
  grid.rect(gp=gpar(col="gray"))
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=j))
      grid.rect()
      grid.text(order[i, j])
      popViewport()
    }
  }
  popViewport()
  grid.text(lab, y=unit(2, "lines"))
  popViewport()
}
pushViewport(viewport(layout=grid.layout(2, 2)))
grid.lshow(1, 1, "(a)", cbind(c(1, 3, 5), c(2, 4, 6)), 3, 2, rep(1, 3), 
  FALSE)
grid.lshow(1, 2, "(b)", cbind(c(6, 4, 2), c(5, 3, 1)), 3, 2, rep(1, 3), 
  FALSE)
grid.lshow(2, 1, "(c)", matrix(c(1, 2), ncol=1), 2, 1, c(2, 1), FALSE)
grid.lshow(2, 2, "(d)", matrix(c(1, 2), ncol=1), 2, 1, c(2, 1), TRUE)
popViewport()



}
figure3.13 <- function() {
grid.lshow <- function(i, j, lab, locs, nrow, ncol, heights, respect) {
  pushViewport(viewport(layout.pos.col=j, layout.pos.row=i))
  pushViewport(viewport(width=unit(1, "npc") - unit(2, "lines"),
               height=unit(1, "npc") - unit(3, "lines"),
	       y=unit(3, "lines"), just="bottom", 
    layout=grid.layout(nrow, ncol, heights=heights, 
      respect=respect)))
  grid.rect(gp=gpar(col="gray"))
  for (i in locs) {
      pushViewport(viewport(layout.pos.row=i$rows, layout.pos.col=i$cols))
      grid.rect()
      grid.text(i$lab)
      popViewport()
  }
  popViewport()
  grid.text(lab, y=unit(2, "lines"))
  popViewport()
}
pushViewport(viewport(layout=grid.layout(2, 2)))
grid.lshow(1, 1, "(a)", 
  list(
    list(rows=1, cols=1, lab=1),
    list(rows=3, cols=1, lab=2)),
  3, 1,
  unit(c(2, 0.5, 1), c("null", "cm", "null")), 
  TRUE)
grid.lshow(1, 2, "(b)", 
  list(
    list(rows=1, cols=1, lab=1),
    list(rows=3, cols=1:2, lab=2),
    list(rows=1, cols=2, lab=3)), 
  3, 2,
  unit(c(2, 0.5, 1), c("null", "cm", "null")), 
  TRUE)
grid.lshow(2, 1, "(c)", 
  list(
    list(rows=1, cols=1, lab=1),
    list(rows=3, cols=1:2, lab=2),
    list(rows=1, cols=2, lab=3)), 
  3, 2,
  unit(c(2, 0.5, 1), c("null", "cm", "null")), 
  cbind(c(0, 0, 1), c(0, 0, 0)))
popViewport()



}
figure3.14 <- function() {
par(mfrow=c(1, 2), mar=c(1, 1, 2, 1))
par(cex=0.7)
x <- 1:10
y <- matrix(sort(rnorm(30)), ncol=3)
plot(x, y[,1], ylim=range(y), ann=FALSE, axes=FALSE, 
     type="l", col="gray")
box(col="gray")

points(x, y[,1])
lines(x, y[,2], col="gray")
points(x, y[,2], pch=2)
lines(x, y[,3], col="gray")
points(x, y[,3], pch=3)


mtext("points() & lines()", side=3, line=0.5)
x <- 1:5
y <- x
plot(x, y, ann=FALSE, axes=FALSE, col="gray", pch=16)
box(col="gray")

text(x[-3], y[-3], c("right", "top", "bottom", "left"), 
     pos=c(4, 3, 1, 2))
text(3, 3, "overlay")

mtext("text()", side=3, line=0.5)



}
figure3.15 <- function() {
t <- seq(60, 360, 30)
x <- cos(t/180*pi)*t/360
y <- sin(t/180*pi)*t/360




source(system.file("extra", "as.raster.R", package="RGraphics"))
rlogo <- pixmap::read.pnm(system.file("pictures/logo.pgm", 
                              package="pixmap")[1])



par(mfrow=c(1, 2), mar=c(1, 1, 2, 1))
par(cex=0.7)

t <- seq(60, 360, 30)
x <- cos(t/180*pi)*t/360
y <- sin(t/180*pi)*t/360


par(mfrow=c(3, 3), mar=rep(1, 4), pty="s")
plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("lines()", side=3, line=.6, cex=.7, family="mono")
lines(x, y)

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("segments()", side=3, line=.6, cex=.7, family="mono")
segments(0, 0, x, y)

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("arrows()", side=3, line=.6, cex=.7, family="mono")
arrows(0, 0, x[-1], y[-1], length=.1)

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("xspline()", side=3, line=.6, cex=.7, family="mono")
xspline(x, y, shape=1)

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("rect()", side=3, line=.6, cex=.7, family="mono")
rect(min(x), min(y), max(x), max(y), col="gray")

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("polygon()", side=3, line=.6, cex=.7, family="mono")
polygon(x, y, col="gray")

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("polypath()", side=3, line=.6, cex=.7, family="mono")
polypath(c(x, NA, .5*x), c(y, NA, .5*y),
         col="gray", rule="evenodd")

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("xspline()", side=3, line=.6, cex=.7, family="mono")
xspline(x, y, shape=1, open=FALSE, col="gray")

plot(x, y, pch=16, col="gray",
     xlim=c(-.6, 1.1), ylim=c(-1.1, .6),
     axes=FALSE, ann=FALSE)
box(col="gray")
mtext("rasterImage()", side=3, line=.6, cex=.7, family="mono")
rasterImage(rlogo,
            x - .07, y - .07,
            x + .07, y + .07,
            interpolate=FALSE)





}
figure3.16 <- function() {
par(mfrow=c(1, 2), mar=c(1, 1, 2, 1), pty="s")
par(cex=0.7)
x <- runif(20, 1, 10)
y <- x + rnorm(20)
plot(x, y, ann=FALSE, axes=FALSE, col="gray", pch=16)
box(col="gray")

lmfit <- lm(y ~ x)
abline(lmfit)
arrows(5, 8, 7, predict(lmfit, data.frame(x=7)),
       length=0.1)
text(5, 8, "Line of best fit", pos=2)

mtext("abline() & arrows()", side=3, line=0.5)
y <- rnorm(50)
hist(y, main="", xlab="", ylab="", axes=FALSE, 
     border="gray", col="light gray")
box(col="gray")
rug(y, ticksize=0.02)

mtext("rug()", side=3, line=0.5)



}
figure3.17 <- function() {
angle <- seq(0, 2*pi, length=13)[-13]
x <- 0.15*cos(angle)
y <- 0.5 + 0.3*sin(angle)
par(mar=rep(0, 4))
plot.new()
box("outer", col="gray")
polygon(0.25 + x, y, col="gray")
text(0.75 + x[c(1, 5, 9)], y[c(1, 5, 9)], "NA", col="gray")
x[c(1, 5, 9)] <- NA
y[c(1, 5, 9)] <- NA
polygon(0.75 + x, y, col="gray")




}
figure3.18 <- function() {
par(mar=c(2, 1, 1, 1))
y1 <- rnorm(100)
y2 <- rnorm(100)

par(mfrow=c(2, 1), xpd=NA)

plot(y1, type="l", axes=FALSE,
     xlab="", ylab="", main="")
box(col="gray")
mtext("Left end of margin", adj=0, side=3)
lines(x=c(20, 20, 40, 40), y=c(-7, max(y1), max(y1), -7), 
      lwd=3, col="gray")

plot(y2, type="l", axes=FALSE,
     xlab="", ylab="", main="")
box(col="gray")
mtext("Right end of margin", adj=1, side=3)
mtext("Label below x=30", at=30, side=1)
lines(x=c(20, 20, 40, 40), y=c(7, min(y2), min(y2), 7), 
      lwd=3, col="gray")




}
figure3.19 <- function() {
par(mfrow=c(2, 1), mar=c(5, 3, 2, 1), cex=0.5, pty="s")
with(iris,
     plot(Sepal.Length, Sepal.Width, 
          pch=as.numeric(Species), cex=1.2))
legend(6.1, 4.4, c("setosa", "versicolor", "virginica"), 
       cex=1.5, pch=1:3)

barplot(VADeaths[1:2,], angle=c(45, 135), density=20, 
        col="gray", names=c("RM", "RF", "UM", "UF"))
legend(0.4, 38, c("55-59", "50-54"), cex=1.5,
       angle=c(135, 45), density=20, fill="gray")




}
figure3.20 <- function() {
par(cex=0.8)
x <- 1:2
y <- runif(2, 0, 100)
par(mar=c(4, 4, 2, 4))
plot(x, y, type="n", xlim=c(0.5, 2.5), ylim=c(-10, 110),
     axes=FALSE, ann=FALSE)

axis(2, at=seq(0, 100, 20))
mtext("Temperature (Centigrade)", side=2, line=3)

axis(1, at=1:2, labels=c("Treatment 1", "Treatment 2"))
axis(4, at=seq(0, 100, 20), labels=seq(0, 100, 20)*9/5 + 32)
mtext("Temperature (Fahrenheit)", side=4, line=3)
box()

segments(x, 0, x, 100, lwd=20)
segments(x, 0, x, 100, lwd=16, col="white")
segments(x, 0, x, y, lwd=16, col="gray")




}
figure3.21 <- function() {
par(mar=rep(1, 4))
plot(0:1, 0:1, type="n", axes=FALSE, ann=FALSE)
usr <- par("usr")
pin <- par("pin")
xcm <- diff(usr[1:2])/(pin[1]*2.54)
ycm <- diff(usr[3:4])/(pin[2]*2.54)

par(xpd=NA)
rect(0 + 0.2*xcm, 0 - 0.2*ycm,
     1 + 0.2*xcm, 1 - 0.2*ycm,
     col="gray", border=NA)

rect(0, 0, 1, 1, col="white")
segments(seq(1, 8, 0.1)*xcm, 0,
         seq(1, 8, 0.1)*xcm, 
         c(rep(c(0.5, rep(0.25, 4), 
                 0.35, rep(0.25, 4)),
               7), 0.5)*ycm)
text(1:8*xcm, 0.6*ycm, 0:7, adj=c(0.5, 0))
text(8.2*xcm, 0.6*ycm, "cm", adj=c(0, 0))




}
figure3.22 <- function() {
layout(matrix(1:2, ncol=1), heights=1:2/6.5)
par(cex=0.7)
drunkenness <- ts(c(3875, 4846, 5128, 5773, 7327, 
                    6688, 5582, 3473, 3186,
                    rep(NA, 51)),
                  start=1912, end=1971)

# Have to copy-and-paste to shrink the mtext text (arggh!)
par(mar=c(5, 6, 2, 4))
plot(drunkenness, lwd=3, col="gray", ann=FALSE, las=2)
mtext("Drunkenness\nRelated Arrests", side=2, line=3.5, cex=0.7)
par(new=TRUE)
plot(nhtemp, ann=FALSE, axes=FALSE)
mtext("Temperature (F)", side=4, line=3, cex=0.7)
title("Using par(new=TRUE) or par(usr=...)")
axis(4)

par(mar=c(5, 4, 4, 2))
with(trees, 
     {
       plot(Height, Volume, pch=3,
            xlab="Height (ft)", 
            ylab=expression(paste("Volume ", (ft^3))))
       symbols(Height, Volume, circles=Girth/12, 
               fg="gray", inches=FALSE, add=TRUE)
     })

mtext("symbols(..., add=TRUE)", font=2, side=3, line=1)



}
figure3.23 <- function() {
xx <- c(1:50)
yy <- rnorm(50)
n <- 50
hline <- 0



xx <- c(1:50)
yy <- rnorm(50)
n <- 50
hline <- 0

par(mfrow=c(2,2), mar=c(3, 3, 1, 1))
plot (yy ~ xx, type="n", axes=FALSE, ann=FALSE)
polygon(c(xx[1], xx, xx[n]), c(min(yy), yy, min(yy)), 
        col="gray", border=NA)

box(col="gray")
plot (yy ~ xx, type="n", axes=FALSE, ann=FALSE)
polygon(c(xx[1], xx, xx[n]), c(min(yy), yy, min(yy)), 
        col="gray", border=NA)

usr <- par("usr")
rect(usr[1], usr[3], usr[2], hline, col="white", border=NA)

box(col="gray")
plot (yy ~ xx, type="n", axes=FALSE, ann=FALSE)
polygon(c(xx[1], xx, xx[n]), c(min(yy), yy, min(yy)), 
        col="gray", border=NA)

usr <- par("usr")
rect(usr[1], usr[3], usr[2], hline, col="white", border=NA)

lines(xx, yy)

box(col="gray")
plot (yy ~ xx, type="n", axes=FALSE, ann=FALSE)
polygon(c(xx[1], xx, xx[n]), c(min(yy), yy, min(yy)), 
        col="gray", border=NA)

usr <- par("usr")
rect(usr[1], usr[3], usr[2], hline, col="white", border=NA)

lines(xx, yy)

abline (h=hline,col="gray")
box()
axis(1)
axis(2) 




}
figure3.24 <- function() {
par(mfrow=c(1, 2), mar=c(3, 3, 1, 1), cex=0.7)
y <- sample(1:10)
midpts <- barplot(y, col=" light gray")
width <- diff(midpts[1:2])/4
left <- rep(midpts, y - 1) - width
right <- rep(midpts, y - 1) + width
heights <- unlist(apply(matrix(y, ncol=10), 
                        2, seq))[-cumsum(y)]
segments(left, heights, right, heights,
         col="white")

with(ToothGrowth, 
     {
       boxplot(len ~ supp, border="gray", 
               col="light gray", boxwex=0.5)
       points(jitter(rep(1:2, each=30), 0.5), 
              unlist(split(len, supp)),
              cex=0.5, pch=16)
     })




}
figure3.25 <- function() {
par(cex=.7)
pairs(iris[1:2], 
      diag.panel=function(x, ...) { 
          boxplot(x, add=TRUE, axes=FALSE,
                  at=mean(par("usr")[1:2])) 
      }, 
      text.panel=function(x, y, labels, ...) { 
          mtext(labels, side=3, line=0) 
      })





}
figure3.26 <- function() {
par(mar=rep(0, 4))
z <- 2 * volcano        
x <- 10 * (1:nrow(z))   
y <- 10 * (1:ncol(z))   
trans <- persp(x, y, z, zlim=c(0, max(z)),
               theta = 150, phi = 12, lwd=.5,
               scale = FALSE, axes=FALSE)

clines <- contourLines(x, y, z)
lapply(clines,
       function(contour) {
           lines(trans3d(contour$x, contour$y, 0, trans))
       })




}
figure3.27 <- function() {
groups <- dimnames(Titanic)[[1]]
males <- Titanic[, 1, 2, 2]
females <- Titanic[, 2, 2, 2]

par(mar=c(0.5, 3, 0.5, 1))

plot.new()
plot.window(xlim=c(-200, 200), ylim=c(-1.5, 4.5))

ticks <- seq(-200, 200, 100)
y <- 1:4
h <- 0.2

lines(rep(0, 2), c(-1.5, 4.5), col="gray")
segments(-200, y, 200, y, lty="dotted")
rect(-males, y-h, 0, y+h, col="dark gray")
rect(0, y-h, females, y+h, col="light gray")
mtext(groups, at=y, adj=1, side=2, las=2)
par(cex.axis=0.5, mex=0.5)
axis(1, at=ticks, labels=abs(ticks), pos=0)

tw <- 1.5*strwidth("females")
rect(-tw, -1-h, 0, -1+h, col="dark gray")
rect(0, -1-h, tw, -1+h, col="light gray")
text(0, -1, "males", pos=2)
text(0, -1, "females", pos=4)

box("inner", col="gray")



}
plot.newclass <- 
  function(x, y=NULL, 
           main="", sub="",
           xlim=NULL, ylim=NULL,
           axes=TRUE, ann=par("ann"),
           col=par("col"),
           ...) {
  xy <- xy.coords(x, y)
  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  plot.new()
  plot.window(xlim, ylim, ...)
  points(xy$x, xy$y, col=col, ...)
  if (axes) {
    axis(1)
    axis(2)
    box()
  }
  if (ann) 
    title(main=main, sub=sub, 
          xlab=xy$xlab, ylab=xy$ylab, ...)
}


figure6.1 <- function() {
pushViewport(viewport(layout=grid.layout(2, 2), gp=gpar(cex=0.6, fill=NA)))
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
pushViewport(plotViewport(c(5, 4, 2, 2)))
pushViewport(dataViewport(pressure$temperature, 
                          pressure$pressure,
                          name="plotRegion"))

grid.points(pressure$temperature, pressure$pressure, 
  gp=gpar(cex=0.5))
grid.rect()
grid.xaxis()
grid.yaxis()
grid.text("temperature", y=unit(-3, "line"))
grid.text("pressure", x=unit(-3, "line"), rot=90)

popViewport(3)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
pushViewport(plotViewport(c(5, 4, 2, 2)))
pushViewport(dataViewport(pressure$temperature, 
                          pressure$pressure,
                          name="plotRegion"))

grid.points(pressure$temperature, pressure$pressure, pch=2, 
  gp=gpar(cex=0.5))
grid.rect()
grid.xaxis()
grid.yaxis()
grid.text("temperature", y=unit(-3, "line"))
grid.text("pressure", x=unit(-3, "line"), rot=90)

popViewport(3)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
pushViewport(plotViewport(c(5, 4, 2, 2)))
pushViewport(dataViewport(pressure$temperature, 
                          pressure$pressure,
                          name="plotRegion"))

grid.points(pressure$temperature, pressure$pressure, pch=2, 
  gp=gpar(cex=0.5))
grid.rect()
grid.xaxis()
grid.yaxis()
grid.text("temperature", y=unit(-3, "line"))
grid.text("pressure", x=unit(-3, "line"), rot=90)

upViewport(2)
grid.rect(gp=gpar(lty="dashed"))

downViewport("plotRegion")
grid.text("Pressure (mm Hg)\nversus\nTemperature (Celsius)",
          x=unit(150, "native"), y=unit(600, "native"))




}
figure6.2 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.circle(x=seq(0.1, 0.9, length=100), 
            y=0.5 + 0.4*sin(seq(0, 2*pi, length=100)),
            r=abs(0.1*cos(seq(0, 2*pi, length=100))))




}
figure6.3 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.circle(c(.1, .3, .4, .6, .7, .9), 
            c(.25, .75), r=unit(1, "mm"),
            gp=gpar(col=NA, fill="gray"))
grid.curve(x1=.1, y1=.25, x2=.3, y2=.75)
grid.curve(x1=.4, y1=.25, x2=.6, y2=.75,
           square=FALSE, ncp=8, curvature=.5)
grid.curve(x1=.7, y1=.25, x2=.9, y2=.75,
           square=FALSE, angle=45, shape=-1)




}
figure6.4 <- function() {
grid.rect(gp=gpar(col="gray"))
angle <- seq(0, 2*pi, length=50)
grid.lines(x=seq(0.1, 0.5, length=50), 
           y=0.5 + 0.3*sin(angle), arrow=arrow())
grid.segments(6:8/10, 0.2, 7:9/10, 0.8,
              arrow=arrow(angle=15, type="closed"))




}
figure6.5 <- function() {
grid.rect(gp=gpar(col="gray"))
angle <- seq(0, 2*pi, length=10)[-10]
grid.polygon(x=0.25 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             gp=gpar(fill="gray"))
grid.polygon(x=0.75 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             id=rep(1:3, each=3),
             gp=gpar(fill="gray"))




}
figure6.6 <- function() {
grid.rect(gp=gpar(col="gray"))
angle <- seq(0, 2*pi, length=10)[-10]
grid.path(x=0.25 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
          gp=gpar(fill="gray"))
grid.path(x=c(0.75 + 0.15*cos(angle), .7, .7, .8, .8),
          y=c(0.5 + 0.3*sin(angle),  .4, .6, .6, .4), 
          id=rep(1:2, c(9, 4)),
          gp=gpar(fill="gray"))




}
figure6.7 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(col="gray")))
grid.text("very snug", 0.4, unit(1, "in"), just=c("left", "bottom"))
grid.lines(x=0.4, y=unit(0:1, "in"), arrow=arrow(ends="both",
                                       length=unit(1, "mm")))
grid.text("1 inch", unit(0.4, "npc") + unit(0.5, "line"), 
  unit(0.5, "in"), rot=90)
grid.lines(x=c(0, 0.4), y=unit(1, "in"), arrow=arrow(ends="both",
                                           length=unit(1, "mm")))
grid.text(unit(0.4, "npc"), 0.2, unit(1, "in") + unit(0.5, "line"))
popViewport()
pushViewport(viewport(gp=gpar(fill=NA)))
grid.rect(x=unit(0.4, "npc"), y=unit(1, "in"),
          width=stringWidth("very snug"), 
          height=unit(1, "line"),
          just=c("left", "bottom"))




}
figure6.8 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(fontsize=10)))
grid.rect(x=0.33, height=0.7, width=0.2, gp=gpar(fill="black"))
grid.rect(x=0.66, height=0.7, width=0.2)
grid.text("grid.rect()", x=0.66, rot=90)
grid.text("grid.rect(gp=gpar(fill=\"black\"))", x=0.33, rot=90, 
  gp=gpar(fontsize=8, col="white"))
popViewport()



}
figure6.9 <- function() {
grid.rect(gp=gpar(col="gray"))
levels <- round(seq(90, 10, length=25))
grays <- paste("gray", c(levels, rev(levels)), sep="")
grid.circle(x=seq(0.1, 0.9, length=100), 
            y=0.5 + 0.4*sin(seq(0, 2*pi, length=100)),
            r=abs(0.1*cos(seq(0, 2*pi, length=100))),
            gp=gpar(col=grays))




}
figure6.10 <- function() {
grid.rect(gp=gpar(col="gray"))
angle <- seq(0, 2*pi, length=11)[-11]
grid.polygon(x=0.25 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             id=rep(1:2, c(7, 3)),
             gp=gpar(fill=c("gray", "white")))
angle[4] <- NA
grid.polygon(x=0.75 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             id=rep(1:2, c(7, 3)),
             gp=gpar(fill=c("gray", "white")))

angle <- seq(0, 2*pi, length=11)[4]
grid.text("NA", x=0.75 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle),
          gp=gpar(col="gray"))



}
figure6.11 <- function() {
vp1 <- 
viewport(x=unit(0.4, "npc"), y=unit(1, "cm"),
         width=stringWidth("very very snug indeed"), 
         height=unit(6, "line"),
         just=c("left", "bottom"))

grid.show.viewport(scale.col="gray", border.fill="white", vp.col="black", vp.fill="gray", vp1)
grid.rect(gp=gpar(col="white", fill=NA, lwd=3))
pushViewport(viewport(.5, .5, .8, .8))
pushViewport(vp1)
grid.rect(gp=gpar(fill=NA))
grid.text("very very snug indeed", 
          gp=gpar(col="white"))
popViewport(2)



}
figure6.12 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))
pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp1"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))




}
figure6.13 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))
pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp1"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))

pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp2"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))




}
figure6.14 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))
pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp1"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))

pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp2"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))

popViewport()
grid.text("bottom-right corner", 
          x=unit(1, "npc") - unit(1, "mm"),
          y=unit(1, "mm"), just=c("right", "bottom"))




}
figure6.15 <- function() {
pushViewport(viewport(gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray"))
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))
pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp1"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))

pushViewport(viewport(width=0.8, height=0.5, angle=10, 
             name="vp2"))
grid.rect()
grid.text("top-left corner", x=unit(1, "mm"),
          y=unit(1, "npc") - unit(1, "mm"), 
          just=c("left", "top"))

popViewport()
grid.text("bottom-right corner", 
          x=unit(1, "npc") - unit(1, "mm"),
          y=unit(1, "mm"), just=c("right", "bottom"))

upViewport()
grid.text("bottom-right corner", 
          x=unit(1, "npc") - unit(1, "mm"),
          y=unit(1, "mm"), just=c("right", "bottom"))
downViewport("vp1")
grid.rect(width=unit(1, "npc") + unit(2, "mm"),
          height=unit(1, "npc") + unit(2, "mm"))




}
figure6.16 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(fill=NA)))
pushViewport(viewport(width=.5, height=.5, clip="on"))
grid.rect()
grid.circle(r=.7, gp=gpar(lwd=20))

pushViewport(viewport(clip="inherit"))
grid.circle(r=.7, gp=gpar(lwd=10, col="gray"))

pushViewport(viewport(clip="off"))
grid.circle(r=.7)
popViewport(3)




}
figure6.17 <- function() {
circText <- function(lab, x, y, suffix) {
    grid.circle(x, y, r=unit(3, "mm"), 
                name=paste(lab, suffix, sep="-"))
    grid.text(lab, x, y,
              gp=if (lab == "ROOT") gpar(cex=.7) else NULL)
}
edge <- function(a, b, angle) {
    grid.segments(grobX(a, angle), grobY(a, angle),
                  grobX(b, 180 + angle), grobY(b, 180 + angle),
                  arrow=arrow(length=unit(2, "mm"),
                    type="closed"),
                  gp=gpar(fill="black"))
}
grid.newpage()
pushViewport(viewport(width=.9, height=.9,
                      layout=grid.layout(2, 2),
                      gp=gpar(cex=.5)))
pushViewport(viewport(layout.pos.col=1, 
                      layout.pos.row=1))
grid.rect(width=.9, height=.9, gp=gpar(col="gray"))
circText("ROOT", .5, .8, 1)
circText("A", .3, .6, 1)
circText("B", .5, .6, 1)
circText("C", .7, .6, 1)
edge("ROOT-1", "A-1", 225)
edge("ROOT-1", "B-1", 270)
edge("ROOT-1", "C-1", 315)
popViewport()
pushViewport(viewport(layout.pos.col=2, 
                      layout.pos.row=1))
grid.rect(width=.9, height=.9, gp=gpar(col="gray"))
circText("ROOT", .5, .8, 2)
circText("A", .5, .6, 2)
circText("B", .5, .4, 2)
circText("C", .5, .2, 2)
edge("ROOT-2", "A-2", 270)
edge("A-2", "B-2", 270)
edge("B-2", "C-2", 270)
popViewport()
pushViewport(viewport(layout.pos.col=1, 
                      layout.pos.row=2))
grid.rect(width=.9, height=.9, gp=gpar(col="gray"))
circText("ROOT", .5, .8, 3)
circText("A", .5, .6, 3)
circText("B", .4, .4, 3)
circText("C", .6, .4, 3)
edge("ROOT-3", "A-3", 270)
edge("A-3", "B-3", 244)
edge("A-3", "C-3", 296)
popViewport()



}
figure6.18 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(fill="gray", fontsize=10)))
grid.text("viewport(gp=gpar(fill=\"gray\"))", y=0.925)
grid.rect(x=0.33, height=0.7, width=0.2)
grid.text("grid.rect()", x=0.33, rot=90)
grid.rect(x=0.66, height=0.7, width=0.2, gp=gpar(fill="black"))
grid.text("grid.rect(gp=gpar(fill=\"black\"))", x=0.66, rot=90, 
  gp=gpar(fontsize=8, col="white"))
popViewport()



}
figure6.19 <- function() {
labelvp <- function(name, col="gray", tcol="white", clipOff=TRUE) {
  seekViewport(name)
  if (clipOff)
    pushViewport(viewport(clip="off"))
  grid.rect(gp=gpar(col=col, lwd=5, fill=NA))
  grid.rect(x=0, y=1, width=unit(1, "strwidth", name) + unit(2, "mm"),
    height=unit(1, "line"), just=c("left", "top"),
    gp=gpar(fill=col, col=NA))
  grid.text(name, x=unit(1, "mm"), y=unit(1, "npc") - unit(1, "mm"),
    just=c("left", "top"), gp=gpar(col=tcol))
  upViewport(0)
}



vplay <- grid.layout(3, 3, 
                     respect=rbind(c(0, 0, 0), 
                                   c(0, 1, 0), 
                                   c(0, 0, 0)))



pushViewport(viewport(width=0.95, height=0.95))
grid.rect(gp=gpar(col="light gray"))
pushViewport(viewport(layout=vplay))

pushViewport(viewport(layout.pos.col=2, name="col2"))
upViewport()
pushViewport(viewport(layout.pos.row=2, name="row2"))

labelvp("col2", "black")
labelvp("row2")



}
figure6.20 <- function() {
unitlay <- 
  grid.layout(3, 3, 
              widths=unit(c(1, 1, 2), 
                          c("in", "null", "null")), 
              heights=unit(c(3, 1, 1), 
                           c("line", "null", "null")))



pushViewport(viewport(gp=gpar(cex=0.8)))
grid.show.layout(unitlay, bg="white", 
                 cell.border="black", cell.fill="gray90", 
                 label.col="black", unit.col="black",
                 newpage=FALSE)
grid.rect(gp=gpar(col="white", lwd=3, fill=NA))
popViewport()



}
figure6.21 <- function() {
gridfun <- function() {
  pushViewport(viewport(layout=grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col=1))
  grid.rect()
  grid.text("black")
  grid.text("&", x=1)
  popViewport()
  pushViewport(viewport(layout.pos.col=2, clip="on"))
  grid.rect(gp=gpar(fill="black"))
  grid.text("white", gp=gpar(col="white"))
  grid.text("&", x=0, gp=gpar(col="white"))
  popViewport(2)
}



grid.rect(gp=gpar(col="gray"))
w <- unit(1, "npc") - unit(15, "mm")
x <- unit.c(unit(5, "mm"),
            unit(5, "mm") + 1/3*w,
            unit(5, "mm") + 1/3*w + unit(5, "mm"),
	    unit(1, "npc") - unit(5, "mm"))
y <- unit.c(unit(5, "mm"),
            unit(5, "mm") + 2/3*w,
            unit(5, "mm") + 2/3*w + unit(5, "mm"),
	    unit(1, "npc") - unit(5, "mm"))
grid.segments(x, 0, x, 1,
  gp=gpar(col="gray", lty="dashed"))
grid.segments(0, y, 1, y,
  gp=gpar(col="gray", lty="dashed"))
pushViewport(
  viewport(
    layout=grid.layout(5, 5, 
                       widths=unit(c(5, 1, 5, 2, 5),
                                   c("mm", "null", "mm",
                                     "null", "mm")),  
                       heights=unit(c(5, 1, 5, 2, 5),
                                    c("mm", "null", "mm",
                                      "null", "mm")))))
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
gridfun()
popViewport()

pushViewport(viewport(layout.pos.col=4, layout.pos.row=4))
gridfun()
popViewport(2)




}
figure6.22 <- function() {

n <- 7
primtest2 <- function(nas, na) {
  angle <- seq(0, 2*pi, length=n+1)[-(n+1)]
  y <- 0.5 + 0.4*sin(angle)
  x <- 0.5 + 0.4*cos(angle)
  if (any(nas))
    grid.text(paste("NA", (1:n)[nas], sep=""),
              x[nas], y[nas], gp=gpar(col="gray"))
  x[nas] <- na
  y[nas] <- na
  grid.polygon(x, y, gp=gpar(fill="light gray", col=NA))
  grid.lines(x, y, arrow=arrow(), gp=gpar(lwd=5))
  grid.move.to(x[1], y[1])
  for (i in 2:n) {
    grid.line.to(x[i], y[i], gp=gpar(col="white"))
  }
}
celltest <- function(r, c, nas, na) {
  pushViewport(viewport(layout.pos.col=c,
                        layout.pos.row=r))
  primtest2(nas, na)
  grid.rect(width=0.9, height=0.9, gp=gpar(col="gray", fill=NA))
  popViewport()
}
cellnas <- function(i) {
  temp <- rep(FALSE, n)
  temp[i] <- TRUE
  temp[n-3+i] <- TRUE
  temp
}
pushViewport(viewport(width=.8, height=.8, 
                      layout=grid.layout(2, 2),
                      gp=gpar(cex=0.7)))
celltest(1, 1, rep(FALSE, n), NA)
celltest(1, 2, cellnas(1), NA)
celltest(2, 1, cellnas(2), NA)
celltest(2, 2, cellnas(3), NA)
popViewport()



}
figure6.23 <- function() {
trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0, ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars,
       panel=function(subscripts, ...) {
           grid.text(paste("n =", length(subscripts)),
                     unit(1, "npc") - unit(1, "mm"),
                     unit(1, "npc") - unit(1, "mm"),
                     just=c("right", "top"))
           panel.xyplot(subscripts=subscripts, ...)
       })

)



}
figure6.24 <- function() {
trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
grid.newpage()
pushViewport(viewport(x=0, width=.4, just="left"))
print(barchart(table(mtcars$gear)),
      newpage=FALSE)
popViewport()
pushViewport(viewport(x=.4, width=.6, just="left"))
print(xyplot(mpg ~ disp, data=mtcars,
             group=gear, 
             auto.key=list(space="right"),
             par.settings=list(superpose.symbol=list(pch=c(1, 3, 16),
                                 fill="white"))),
      newpage=FALSE)
popViewport()



}
figure6.25 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



print(
ggplot(mtcars2, aes(x=disp, y=mpg)) +
    geom_point()

)
downViewport("panel.3-4-3-4")
grid.text(paste("n =", nrow(mtcars2)),
          x=unit(1, "npc") - unit(1, "mm"), 
          y=unit(1, "npc") - unit(1, "mm"),
          just=c("right", "top"))




}
figure6.26 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



grid.newpage()
pushViewport(viewport(x=0, width=1/3, just="left"))
print(ggplot(mtcars2, aes(x=trans)) + 
      geom_bar(),
      newpage=FALSE)
popViewport()
pushViewport(viewport(x=1/3, width=2/3, just="left"))
print(ggplot(mtcars2, aes(x=disp, y=mpg)) +
      geom_point(aes(color=trans)) +
      scale_color_manual(values=gray(2:1/3)),
      newpage=FALSE)
popViewport()



}
figure1.2 <- function() {


#
#  Comment:
# 
#  Examples of the use of standard high-level plotting functions.
# 
#  In each case, extra output is also added using low-level 
#  plotting functions.
#


par(mfrow=c(3, 2))

# Scatterplot
x <- c(0.5, 2, 4, 8, 12, 16)
y1 <- c(1, 1.3, 1.9, 3.4, 3.9, 4.8)
y2 <- c(4, .8, .5, .45, .4, .3)
par(las=1, mar=c(4, 4, 2, 4), cex=.7)
plot.new()
plot.window(range(x), c(0, 6))
lines(x, y1)
lines(x, y2)
points(x, y1, pch=16, cex=2)
points(x, y2, pch=21, bg="white", cex=2)
par(col="gray50", fg="gray50", col.axis="gray50")
axis(1, at=seq(0, 16, 4))
axis(2, at=seq(0, 6, 2))
axis(4, at=seq(0, 6, 2))
box(bty="u")
mtext("Travel Time (s)", side=1, line=2, cex=0.8)
mtext("Responses per Travel", side=2, line=2, las=0, cex=0.8)
mtext("Responses per Second", side=4, line=2, las=0, cex=0.8)
text(4, 5, "Bird 131")
par(mar=c(5.1, 4.1, 4.1, 2.1), col="black", fg="black", col.axis="black")

# Histogram
# Random data
Y <- rnorm(50)
# Make sure no Y exceed [-3.5, 3.5]
Y[Y < -3.5 | Y > 3.5] <- NA
x <- seq(-3.5, 3.5, .1)
dn <- dnorm(x)
par(mar=c(4.5, 4.1, 3.1, 0))
hist(Y, breaks=seq(-3.5, 3.5), ylim=c(0, 0.5), 
     col="gray80", freq=FALSE)
lines(x, dnorm(x), lwd=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Barplot
# Modified from example(barplot)
par(mar=c(2, 3.1, 2, 2.1))
midpts <- barplot(VADeaths, 
                  col=gray(0.5 + 1:5/12), 
                  names=rep("", 4))
mtext(sub(" ", "\n", colnames(VADeaths)),
      at=midpts, side=1, line=0.5, cex=0.5)
text(rep(midpts, each=5), apply(VADeaths, 2, cumsum) - VADeaths/2,
     VADeaths, 
     col=rep(c("white", "black"), times=2:3), 
     cex=0.8)
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Boxplot
# Modified example(boxplot) - itself from suggestion by Roger Bivand
par(mar=c(3, 4.1, 2, 0))
     boxplot(len ~ dose, data = ToothGrowth,
             boxwex = 0.25, at = 1:3 - 0.2,
             subset= supp == "VC", col="gray90",
             xlab="",
             ylab="tooth length", ylim=c(0,35))
     mtext("Vitamin C dose (mg)", side=1, line=2.5, cex=0.8)
     boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
             boxwex = 0.25, at = 1:3 + 0.2,

             subset= supp == "OJ")
     legend(1.5, 9, c("Ascorbic acid", "Orange juice"), 
            fill = c("gray90", "gray70"), 
            bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Persp
# Almost exactly example(persp)
    x <- seq(-10, 10, length= 30)
     y <- x
     f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
     z <- outer(x, y, f)
     z[is.na(z)] <- 1
# 0.5 to include z axis label
par(mar=c(0, 0.5, 0, 0), lwd=0.5)
     persp(x, y, z, theta = 30, phi = 30, 
 
           expand = 0.5)
par(mar=c(5.1, 4.1, 4.1, 2.1), lwd=1)

# Piechart
# Example 4 from help(pie)
par(mar=c(0, 2, 1, 2), xpd=FALSE, cex=0.5)
     pie.sales <- c(0.12, 0.3, 0.26, 0.16, 0.04, 0.12)
     names(pie.sales) <- c("Blueberry", "Cherry",
         "Apple", "Boston Cream", "Other", "Vanilla")
     pie(pie.sales, col = gray(seq(0.4,1.0,length=6))) 




}
figure1.3 <- function() {

#
# Comment:
#
# A sophisticated example of adding further output to a basic plot.
# 
# Most of the functions defined are just for calculating values
# relevant to the data analysis.  
# 
# The function plotPars() is the one of interest for seeing how
# the drawing of the plot is done.
#


params <- function(N, breaks, p=seq(0.001, 1, length=100)) {
  list(N=N, T=1/breaks, p=p, q=1-p)
}

pdfcomp <- function(comp, params) {
  n <- params$T
  p <- params$p
  q <- params$q
  y <- round(comp/n)
  choose(n, comp)*p^comp*q^(n-comp) / (1 - q^n)
}

# Expected num sherds (for a vessel) [=completeness]
expcomp <- function(params) {
  params$T*params$p/(1-params$q^params$T)
}

# Variance of num sherds (for a vessel)
varcomp <- function(params) {
  n <- params$T
  p <- params$p
  q <- params$q
  # From Johnson & Kotz
  (n*p*q / (1 - q^n)) - (n^2*p^2*q^n / (1 - q^n)^2)
  # n^2 times Thomas Yee's formula
  # n^2*((p*(1 + p*(n - 1)) / (n*(1 - q^n))) - (p^2 / (1 - q^n)^2))
}

# Expected value of completeness (for a sample of vessels)
expmeancomp <- function(params) {
  expcomp(params)
}

# Variance of completeness (for a sample of vessels)
# Use the expected number of vessels in sample as denominator
varmeancomp <- function(params) {
  varcomp(params)/(numvess(params))
}

numvess <- function(params) {
  params$N*(1-params$q^params$T)
}

ecomp <- function(p, T, comp) {
  q <- 1 - p
  T*p/(1 - q^T) - comp
}

estN <- function(comp, broke, n) {
  T <- 1/broke
  n / (1 - (1 - uniroot(ecomp, c(0.00001, 1), T=T, comp=comp)$root)^T)
}

nvessscale <- function(params, xlim, ylim, new=TRUE) {
  if (new)
    par(new=TRUE)
  plot(0:1, c(1, params$N), type="n", axes=!new, ann=FALSE,
       xlim=xlim, ylim=ylim)
}

compscale <- function(params, xlim, ylim, new=TRUE) {
  if (new)
    par(new=TRUE)
  plot(0:1, c(1, params$T), type="n", axes=!new, ann=FALSE,
       xlim=xlim, ylim=ylim)
}

lowerCI <- function(p, N, breaks, lb) {
  params <- params(N, breaks, p)
  expmeancomp(params) - 2*sqrt(varmeancomp(params)) - lb
}

upperCI <- function(p, N, breaks, lb) {
  params <- params(N, breaks, p)
  expmeancomp(params) + 2*sqrt(varmeancomp(params)) - lb
}

critP <- function(comp, params) {
  c(uniroot(lowerCI, c(0.00001, 1), N=params$N,
            breaks=1/params$T, lb=max(comp))$root,
    if (upperCI(0.00001, params$N, 1/params$T, min(comp)) > 0) 0
    else uniroot(upperCI, c(0.00001, 1), N=params$N,
                 breaks=1/params$T, lb=min(comp))$root)
}

anncomp <- function(params, comp, xlim, ylim, cylim) {
  cp <- critP(comp, params)
  nv <- numvess(params(params$N, 1/params$T, cp))
  nvessscale(params, xlim, ylim)
  polygon(c(cp[2], cp[2], 0, 0, cp[1], cp[1]),
          c(0, nv[2], nv[2], nv[1], nv[1], 0),
          col="gray90", border=NA)
  text(0, nv[1], paste(round(nv[1]),
                       " (", round(100*nv[1]/params$N), "%)", sep=""),
       adj=c(0, 0), col="gray")
  text(0, nv[2], paste(round(nv[2]), 
                       " (", round(100*nv[2]/params$N), "%)", sep=""),
       adj=c(0, 1), col="gray")
  compscale(params, xlim, cylim)
  segments(1, min(comp), cp[2], comp, col="gray")
  segments(1, max(comp), cp[1], comp, col="gray")
  text(1, comp, paste(comp, collapse="-"), adj=c(1, 0), col="gray")
}

plotPars <- function(params, comp, xlim=NULL, ylim=NULL) {
  mean <- expmeancomp(params)
  var <- 2*sqrt(varmeancomp(params))
  lb <- mean - var
  ub <- mean + var
  par(mar=c(5, 4, 4, 4))
  if (is.null(ylim))
    cylim <- ylim
  else
    cylim <- c(1 + ((ylim[1] - 1)/(params$N - 1))*(params$T - 1),
               1 + ((ylim[2] - 1)/(params$N - 1))*(params$T - 1))
  nvessscale(params, xlim, ylim, new=FALSE)
  compscale(params, xlim, cylim)
  polygon(c(params$p, rev(params$p)), c(lb, rev(ub)),
          col="gray90", border=NA)
  anncomp(params, comp, xlim, ylim, cylim)
  nvessscale(params, xlim, ylim)
  mtext("Number of Vessels", side=2, line=3)
  mtext("Sampling Fraction", side=1, line=3)
  lines(params$p, numvess(params))
  par(new=TRUE)
  compscale(params, xlim, cylim)
  mtext("Completeness", side=4, line=3)
  axis(4)
  lines(params$p, mean, lty="dashed")
  lines(params$p, lb, lty="dotted")
  lines(params$p, ub, lty="dotted")
  mtext(paste("N = ", round(params$N),
              "     brokenness = ", round(1/params$T, 3), sep=""),
        side=3, line=2)
}

par(cex=0.8, mar=c(3, 3, 3, 3))
p6 <- params(estN(1.2, 0.5, 200), 0.5)
plotPars(p6, 1.2)
nvessscale(p6, NULL, NULL)
pcrit <- 1 - (1 - 200/estN(1.2, 0.5, 200))^(1/p6$T)
lines(c(0, pcrit), c(200, 200))
lines(c(pcrit, pcrit), c(200, 0))



}
figure1.4 <- function() {
# Produce a plot of tiger populations with picture as background
# Source: http://www.globaltiger.org/population.htm
year <- c(1993, 1996, 1998, 2001)
minpop <- c(20, 50, 50, 115)
maxpop <- c(50, 240, 240, 150)


tiger <- grImport::readPicture(system.file("extra", "tiger.ps.xml", 
                                 package="RGraphics"))[-1]



source(system.file("extra", "grayify.R", package="RGraphics"))

# grid.newpage()
pushViewport(plotViewport(c(3, 2, 2, 1)),
             viewport(xscale=c(1991, 2003), yscale=c(0, 250)))
grid.rect()
# tiger backdrop in gray
grImport::grid.picture(tiger, x=0.45, FUN=grayify, min=.8)
grid.xaxis(at=year, gp=gpar(cex=0.7))
grid.yaxis(gp=gpar(cex=0.7))
# black bars
grid.rect(x=unit(year, "native"), y=0,
          width=unit(1, "native"), height=unit(maxpop, "native"),
          just="bottom", gp=gpar(fill="black"))
# tiger in bars
tigerGrob <- grImport::pictureGrob(tiger, x=0.45, 
 FUN=grImport::grobify)
# Start from 2 because bar 1 does not overlap with tiger
for (i in 2:length(year)) {
    grid.clip(x=unit(year[i], "native"), y=0,
              width=unit(1, "native"), height=unit(maxpop[i], "native"),
              just="bottom")
    # tiger backdrop (shift slightly to left so get one eye in one bar)
    grid.draw(tigerGrob)
}
grid.clip()
# redo bar borders
grid.rect(x=unit(year, "native"), y=0,
          width=unit(1, "native"), height=unit(maxpop, "native"),
          just="bottom", gp=gpar(fill=NA))
grid.text("Estimated Population (max.) of Bengal Tigers\n(in Bhutan)",
          y=unit(1, "npc") + unit(1, "lines"))
popViewport(2)



}
figure1.5 <- function() {

#
# Comment:
#
# A slightly modified version of Figure 1.1 from 
# Cleveland's book "Visualizing Data"
#



trellis.par.set(list(fontsize=list(text=6),
	             par.xlab.text=list(cex=1.5),
                     add.text=list(cex=1.5),
                     superpose.symbol=list(cex=.5)))
key <- simpleKey(levels(lattice::barley$year), space = "right")
key$text$cex <- 1.5
print(
     dotplot(variety ~ yield | site, data = lattice::barley, groups = year,
             key = key,
             xlab = "Barley Yield (bushels/acre) ",
             aspect=0.5, layout = c(1,6), ylab=NULL)
)



}
figure1.6 <- function() {

#
# Comment:
#
# Inspired by Figure 3.3 from 
# Wickham's book "ggplot2"
#


print(
ggplot(data=ggplot2::mpg, aes(x=displ, y=hwy, shape=factor(cyl))) + 
    geom_point() +
    stat_smooth(method="lm", colour="black") +
    scale_shape_manual(values=c(1, 16, 3, 17)) + 
    theme_bw() 
)



}
figure1.7 <- function() {

#
# Comment:
#
# A bit of mucking around is required to get the second (whole-world)
# map positioned correctly;  this provides an example of calling a 
# plotting function to perform calculations but do no drawing (see the
# second call to the map() function).
#
# Makes use of the "maps", "mapdata", and "mapproj" packages to draw the maps.
#




par(mar=rep(1, 4))
maps::map("nzHires", fill=TRUE, col="gray80",
    regions=c("North Island", "South Island", "Stewart Island"))
points(174.75, -36.87, pch=16, cex=2,
       col=rgb(0,0,0,.5))
arrows(172, -36.87, 174, -36.87, lwd=3)
text(172, -36.87, "Auckland  ", adj=1, cex=2)
# mini world map as guide
maplocs <- maps::map(projection="sp_mercator", wrap=TRUE, lwd=0.1, 
               col="gray", ylim=c(-60, 75),
               interior=FALSE, orientation=c(90, 180, 0), add=TRUE,
               plot=FALSE)
xrange <- range(maplocs$x, na.rm=TRUE)
yrange <- range(maplocs$y, na.rm=TRUE)
aspect <- abs(diff(yrange))/abs(diff(xrange))
# customised to 6.5 by 4.5 figure size
par(fig=c(0.99 - 0.5, 0.99, 0.01, 0.01 + 0.5*aspect*4.5/6.5), 
    mar=rep(0, 4), new=TRUE)
plot.new()
plot.window(xlim=xrange,
            ylim=yrange)
maps::map(projection="sp_mercator", wrap=TRUE, lwd=0.5, ylim=c(-60, 75),
    interior=FALSE, orientation=c(90, 180, 0), add=TRUE)
symbols(-.13, -0.8, circles=1, inches=0.1, add=TRUE)
box()



}
figure1.8 <- function() {

quantmod::getSymbols("YHOO")
quantmod::chartSeries(YHOO, subset='last 4 months'
 )




}
figure1.9 <- function() {


# CLASSIFICATION
# fitting

data("GlaucomaM", envir=environment())
glau <- GlaucomaM
levels(glau$Class) <- c("glau", "norm")
fm.class <- party::ctree(Class ~ ., data = glau)

# visualization
pushViewport(viewport(gp=gpar(cex=0.6)))
plot(fm.class, new=FALSE, terminal.panel=myNode)
popViewport()



}
figure1.10 <- function() {

#
# Comment:
#
# Some simple ideas as a basis for meta-analysis plots.
# 
# The code is modular so that something similar could be achieved
# with different data quite simply.  The actual drawing for these data
# only occurs in the last 10 or so lines of code.
#


# The horizontal gap between columns with content
colgap <- unit(3, "mm")

# The data for column 1
# 
# Of course, many other possible ways to represent the data
# One advantage with this way is that col1$labels can be used
# directly in the calculation of the column widths for the
# main table (see below)
#
# NOTE:  textGrobs are used here so that the fontface (bold in
# some cases) is associated with the label.  In this way, the
# calculation of column widths takes into account the font face.
col1 <- list(labels=
             list(textGrob("Centre", x=0, just="left",
                           gp=gpar(fontface="bold", col="white")),
                  textGrob("Thailand", x=0, just="left"),
                  textGrob("Philippines", x=0, just="left"),
                  textGrob("All in situ", x=0, just="left",
                           gp=gpar(fontface="bold.italic")),
                  textGrob("Colombia", x=0, just="left"),
                  textGrob("Spain", x=0, just="left"),
                  textGrob("All invasive", x=0, just="left",
                           gp=gpar(fontface="bold.italic")),
                  textGrob("All", x=0, just="left",
                           gp=gpar(fontface="bold"))),
             rows=c(1, 5, 6, 8, 11, 12, 14, 16))

# Labels in col 1 which are not used to calculate the
# column width (they spill over into col 2)
col1plus <- list(labels=
                 list(textGrob("Carcinoma in situ", x=0, just="left",
                               gp=gpar(fontface="bold.italic")),
                      textGrob("Invasive cancer", x=0, just="left",
                               gp=gpar(fontface="bold.italic"))),
                 rows=c(4, 10))

# Data for column 2
col2 <- list(labels=
             list(textGrob("Cases", x=1, just="right",
                           gp=gpar(fontface="bold", col="white")),
                  textGrob("327", x=1, just="right"),
                  textGrob("319", x=1, just="right"),
                  textGrob("1462", x=1, just="right",
                           gp=gpar(fontface="bold")),
                  textGrob("96", x=1, just="right"),
                  textGrob("115", x=1, just="right"),
                  textGrob("211", x=1, just="right",
                           gp=gpar(fontface="bold")),
                  textGrob("1673", x=1, just="right",
                           gp=gpar(fontface="bold"))),
             rows=c(1, 5, 6, 8, 11, 12, 14, 16))

# Data for column 3 (width specified as a physical size below)
col3 <- list(OR=c(0.72, 1.27, 1.17, 2.97, 1.86, 2.01, 1.20),
             LL=c(0.52, 0.87, 1.03, 1.42, 0.46, 1.09, 1.07),
             UL=c(1.00, 1.85, 1.32, 6.21, 7.51, 3.71, 1.35),
             rows=c(5, 6, 8, 11, 12, 14, 16),
             # "s" means summary, "n" means normal
             type=c("n", "n", "s", "n", "n", "s", "s"))

# Sizes of boxes
information <- sqrt(1 / ((log(col3$UL) - log(col3$OR))/1.96))
col3$sizes <- information/max(information)

# Width of column 3
col3width <- unit(1.5, "inches")

# Range on the x-axis for column 3
col3$range <- c(0, 4)

# Function to draw a cell in a text column
drawLabelCol <- function(col, j) {
  for (i in 1:length(col$rows)) {
    pushViewport(viewport(layout.pos.row=col$rows[i], layout.pos.col=j))
    # Labels are grobs containing their location so just
    # have to grid.draw() them
    grid.draw(col$labels[[i]])
    popViewport()
  }
}

# Function to draw a non-summary rect-plus-CI
drawNormalCI <- function(LL, OR, UL, size) {
  # NOTE the use of "native" units to position relative to
  # the x-axis scale, and "snpc" units to size relative to
  # the height of the row
  # ("snpc" stands for "square normalised parent coordinates"
  #  which means that the value is calculated as a proportion
  #  of the width and height of the current viewport and the
  #  physically smaller of these is used)
  grid.rect(x=unit(OR, "native"),
            width=unit(size, "snpc"), height=unit(size, "snpc"),
            gp=gpar(fill="black"))
  # Draw arrow if exceed col range
  # convertX() used to convert between coordinate systems
  if (convertX(unit(UL, "native"), "npc", valueOnly=TRUE) > 1)
    grid.lines(x=unit(c(LL, 1), c("native", "npc")), y=.5,
               arrow=arrow(length=unit(0.05, "inches")))
  else {
    # Draw line white if totally inside rect
    lineCol <- if ((convertX(unit(OR, "native") + unit(0.5*size, "lines"),
                             "native", valueOnly=TRUE) > UL) &&
                   (convertX(unit(OR, "native") - unit(0.5*size, "lines"),
                             "native", valueOnly=TRUE) < LL))
      "white"
    else
      "black"
    grid.lines(x=unit(c(LL, UL), "native"), y=0.5,
               gp=gpar(col=lineCol))
  }
}

# Function to draw a summary "diamond"
drawSummaryCI <- function(LL, OR, UL, size) {
  # Not sure how to calc the heights of the diamonds so
  # I'm just using half the height of the equivalent rect
  grid.polygon(x=unit(c(LL, OR, UL, OR), "native"),
               y=unit(0.5 + c(0, 0.25*size, 0, -0.25*size), "npc"))
}

# Function to draw a "data" column
drawDataCol <- function(col, j) {
  pushViewport(viewport(layout.pos.col=j, xscale=col$range))
  grid.lines(x=unit(1, "native"), y=0:1)
  # Assume that last value in col is "All"
  grid.lines(x=unit(col$OR[length(col$OR)], "native"),
             y=0:1, gp=gpar(lty="dashed"))
  grid.xaxis(gp=gpar(cex=0.6))
  grid.text("OR", y=unit(-2, "lines"))
  popViewport()
  for (i in 1:length(col$rows)) {
    pushViewport(viewport(layout.pos.row=col$rows[i], layout.pos.col=j,
                          xscale=col$range))
    if (col$type[i] == "n")
      drawNormalCI(col$LL[i], col$OR[i], col$UL[i], col$sizes[i])
    else
      drawSummaryCI(col$LL[i], col$OR[i], col$UL[i], col$sizes[i])
    popViewport()
  }
}

# Draw the table
#
# The table is just a big layout
#
# All rows are the height of 1 line of text
# 
# Widths of column 1 and 2 are based on widths of labels in
# col$labels and col2$labels 
pushViewport(viewport(layout=grid.layout(16, 5,
                        widths=
                        unit.c(max(unit(rep(1, 8), "grobwidth", col1$labels)),
                               colgap,
                               max(unit(rep(1, 8), "grobwidth", col2$labels)),
                               colgap,
                               col3width),
                        heights=unit(c(1, 0, rep(1, 14)), "lines"))))
pushViewport(viewport(layout.pos.row=1))
grid.rect(gp=gpar(col=NA, fill="black"))
popViewport()
for (i in c(8, 14, 16)) {
    pushViewport(viewport(layout.pos.row=i))
    grid.rect(gp=gpar(col=NA, fill="gray80"))
    popViewport()
}
drawLabelCol(col1, 1)
drawLabelCol(col1plus, 1)
drawLabelCol(col2, 3)
drawDataCol(col3, 5)
popViewport()
                          



}
figure1.11 <- function() {

#
# Comment:
#
# Code by Arden Miller (Department of Statistics, The University of Auckland).
# 
# Lots of coordinate transformations being done "by hand".
# This code is not really reusable;  just a demonstration that very 
# pretty results are possible if you're sufficiently keen.
#


par(mfrow=c(2, 1), pty="s", mar=rep(1, 4)) 
# Create plotting region and plot outer circle
plot(c(-1.1, 1.2), c(-1.1, 1.2),
     type="n", xlab="", ylab="", 
     xaxt="n", yaxt="n", cex.lab=2.5)
angs <- seq(0, 2*pi, length=500)
XX <- sin(angs)
YY <- cos(angs)
lines(XX, YY, type="l")

# Set constants
phi1 <- pi*2/9
k1 <- sin(phi1)
k2 <- cos(phi1)

# Create gray regions
obsphi <- pi/12
lambdas <- seq(-pi, pi, length=500)
xx <- cos(pi/2 - obsphi)*sin(lambdas)
yy <- k2*sin(pi/2 - obsphi)-k1 * cos(pi/2 - obsphi)*cos(lambdas)
polygon(xx, yy, col="gray")
lines(xx, yy, lwd=2)
theta1sA <- seq(-obsphi, obsphi, length=500)
theta2sA <- acos(cos(obsphi)/cos(theta1sA))
theta1sB <- seq(obsphi, -obsphi, length=500)
theta2sB <-  -acos(cos(obsphi)/cos(theta1sB))
theta1s <- c(theta1sA, theta1sB)
theta2s <- c(theta2sA, theta2sB)
xx <- cos(theta1s)*sin(theta2s+pi/4)
yy <- k2*sin(theta1s)-k1*cos(theta1s)*cos(theta2s+pi/4)
polygon(xx, yy, col="gray")
lines(xx, yy, lwd=2)
xx <- cos(theta1s)*sin(theta2s-pi/4)
yy <- k2*sin(theta1s)-k1*cos(theta1s)*cos(theta2s-pi/4)
polygon(xx, yy, col="gray")
lines(xx, yy, lwd=2)

# Plot longitudes
vals <- seq(0, 7, 1)*pi/8
for(lambda in vals){
sl <- sin(lambda)
cl <- cos(lambda)
phi <- atan(((0-1)*k2*cl)/(k1))
angs <- seq(phi, pi+phi, length=500)
xx <- cos(angs)*sl
yy <- k2*sin(angs)-k1*cos(angs)*cl
lines(xx, yy, lwd=.5)
}

# Grey out polar cap
phi <- 5.6*pi/12
lambdas <- seq(-pi, pi, length=500)
xx <- cos(phi)*sin(lambdas)
yy <- k2*sin(phi)-k1 * cos(phi)*cos(lambdas)
polygon(xx, yy, col="gray")

# Plot Latitudes
vals2 <- seq(-2.8, 5.6, 1.4)*pi/12
for(phi in vals2){
  if (k1*sin(phi) > k2 * cos(phi)) 
    crit <- pi 
  else 
    crit <- acos((-k1*sin(phi))/(k2*cos(phi)))
  lambdas <- seq(-crit, crit, length=500)
  xx <- cos(phi)*sin(lambdas)
  yy <- k2*sin(phi)-k1 * cos(phi)*cos(lambdas)
  lines(xx, yy, lwd=.5)
}


# Plots axes and label
lines(c(0.00, 0.00), c(k2*sin(pi/2), 1.11), lwd=4)
lines(c(0.00, 0.00), c(-1, -1.12), lwd=4)
a2x <- sin(-pi/4)
a2y <- cos(-pi/4)*(-k1)
lines(c(a2x, 1.5*a2x), c(a2y, 1.5*a2y), lwd=4)
k <- sqrt(a2x^2+a2y^2)
lines(c(-a2x/k, 1.2*(-a2x/k)), c(-a2y/k, 1.2*(-a2y/k)), lwd=4)
a3x <- sin(pi/4)
a3y <- cos(pi/4)*(-k1)
lines(c(a3x, 1.5*a3x), c(a3y, 1.5*a3y), lwd=4)
k <- sqrt(a3x^2+a3y^2)
lines(c(-a3x/k, 1.2*(-a3x/k)), c(-a3y/k, 1.2*(-a3y/k)), lwd=4)
text(0.1, 1.12, expression(bold(X[1])))
text(-1.07, -.85, expression(bold(X[2])))
text(1.11, -.85, expression(bold(X[3])))

# set plot region and draw outer circle
plot(c(-1.1, 1.2),  c(-1.1, 1.2),
     type="n", xlab="", ylab="", 
     xaxt="n", yaxt="n", cex.lab=2.5)
angs <- seq(0, 2*pi, length=500)
XX <- sin(angs)
YY <- cos(angs)
lines(XX, YY, type="l")

# set constants
phi1 <- pi*2/9
k1 <- sin(phi1)
k2 <- cos(phi1)
obsphi <- pi/24

# create X2X3 gray region and plot boundary
crit <- acos((-k1*sin(obsphi))/(k2 * cos(obsphi)))
lambdas <- seq(-crit, crit, length=500)
xx1 <- cos(obsphi)*sin(lambdas)
yy1 <- k2*sin(obsphi)-k1 * cos(obsphi)*cos(lambdas)
obsphi <-  -pi/24
crit <- acos((-k1*sin(obsphi))/(k2 * cos(obsphi)))
lambdas <- seq(crit, -crit, length=500)
xx3 <- cos(obsphi)*sin(lambdas)
yy3 <- k2*sin(obsphi)-k1 * cos(obsphi)*cos(lambdas)
ang1 <-  atan(xx1[500]/yy1[500])
ang2 <- pi+atan(xx3[1]/yy3[1])
angs <- seq(ang1, ang2, length=50)
xx2 <- sin(angs)
yy2 <- cos(angs)
ang4 <-  atan(xx1[1]/yy1[1])
ang3 <-  -pi+ atan(xx3[500]/yy3[500])
angs <- seq(ang3, ang4, length=50)
xx4 <- sin(angs)
yy4 <- cos(angs)
xxA <- c(xx1, xx2, xx3, xx4)
yyA <- c(yy1, yy2, yy3, yy4)
polygon(xxA, yyA, border="gray", col="gray")
xx1A <- xx1
yy1A <- yy1
xx3A <- xx3
yy3A <- yy3

# create X1X3 gray region and plot boundary
obsphi <- pi/24
crit <- pi/2-obsphi
theta1sA <- c(seq(-crit, crit/2, length=200), seq(crit/2, crit, length=500))
theta2sA <- asin(cos(crit)/cos(theta1sA))
theta1sB <- seq(crit, crit/2, length=500)
theta2sB <-  pi-asin(cos(crit)/cos(theta1sB))
theta1s <- c(theta1sA, theta1sB)
theta2s <- c(theta2sA, theta2sB)
vals <- k1*sin(theta1s)+k2*cos(theta1s)*cos(theta2s+pi/4)
xx1 <- cos(theta1s[vals>=0])*sin(theta2s[vals>=0]+pi/4)
yy1 <- k2*sin(theta1s[vals>=0])-k1*cos(theta1s[vals>=0])*cos(theta2s[vals>=0]+pi/4)
theta2s <-  -theta2s
vals <- k1*sin(theta1s)+k2*cos(theta1s)*cos(theta2s+pi/4)
xx3 <- cos(theta1s[vals>=0])*sin(theta2s[vals>=0]+pi/4)
yy3 <- k2*sin(theta1s[vals>=0])-k1*cos(theta1s[vals>=0])*cos(theta2s[vals>=0]+pi/4)
rev <- seq(length(xx3), 1, -1)
xx3 <- xx3[rev]
yy3 <- yy3[rev]
ang1 <-  pi+atan(xx1[length(xx1)]/yy1[length(yy1)])
ang2 <-  pi+atan(xx3[1]/yy3[1])
angs <- seq(ang1, ang2, length=50)
xx2 <- sin(angs)
yy2 <- cos(angs)
ang4 <-  pi+atan(xx1[1]/yy1[1])
ang3 <-  pi+atan(xx3[length(xx3)]/yy3[length(yy3)])
angs <- seq(ang3, ang4, length=50)
xx4 <- sin(angs)
yy4 <- cos(angs)
xxB <- c(xx1, -xx2, xx3, xx4)
yyB <- c(yy1, -yy2, yy3, yy4)
polygon(xxB, yyB, border="gray", col="gray")
xx1B <- xx1
yy1B <- yy1
xx3B <- xx3
yy3B <- yy3

# create X1X2 gray region and plot boundary
vals <- k1*sin(theta1s)+k2*cos(theta1s)*cos(theta2s-pi/4)
xx1 <- cos(theta1s[vals>=0])*sin(theta2s[vals>=0]-pi/4)
yy1 <- k2*sin(theta1s[vals>=0])-k1*cos(theta1s[vals>=0])*cos(theta2s[vals>=0]-pi/4)
theta2s <-  -theta2s
vals <- k1*sin(theta1s)+k2*cos(theta1s)*cos(theta2s-pi/4)
xx3 <- cos(theta1s[vals>=0])*sin(theta2s[vals>=0]-pi/4)
yy3 <- k2*sin(theta1s[vals>=0])-k1*cos(theta1s[vals>=0])*cos(theta2s[vals>=0]-pi/4)
rev <- seq(length(xx3), 1, -1)
xx3 <- xx3[rev]
yy3 <- yy3[rev]
ang1 <-  pi+atan(xx1[length(xx1)]/yy1[length(yy1)])
ang2 <-  pi+atan(xx3[1]/yy3[1])
angs <- seq(ang1, ang2, length=50)
xx2 <- sin(angs)
yy2 <- cos(angs)
ang4 <-  pi+atan(xx1[1]/yy1[1])
ang3 <-  pi+atan(xx3[length(xx3)]/yy3[length(yy3)])
angs <- seq(ang3, ang4, length=50)
xx4 <- sin(angs)
yy4 <- cos(angs)
xx <- c(xx1, -xx2, xx3, xx4)
yy <- c(yy1, -yy2, yy3, yy4)
polygon(xx, yy, border="gray", col="gray")
xx1C <- xx1
yy1C <- yy1
xx3C <- xx3
yy3C <- yy3


# plot boundaries to gray regions
lines(xx1C[2:45], yy1C[2:45], lwd=2)
lines(xx1C[69:583], yy1C[69:583], lwd=2)
lines(xx1C[660:1080], yy1C[660:1080], lwd=2)
lines(xx3C[13:455], yy3C[13:455], lwd=2)
lines(xx3C[538:1055], yy3C[538:1055], lwd=2)
lines(xx3C[1079:1135], yy3C[1079:1135], lwd=2)
lines(xx1A[6:113], yy1A[6:113], lwd=2)
lines(xx1A[153:346], yy1A[153:346], lwd=2)
lines(xx1A[389:484], yy1A[389:484], lwd=2)
lines(xx3A[1:93], yy3A[1:93], lwd=2)
lines(xx3A[140:362], yy3A[140:362], lwd=2)
lines(xx3A[408:497], yy3A[408:497], lwd=2)
lines(xx1B[2:45], yy1B[2:45], lwd=2)
lines(xx1B[69:583], yy1B[69:583], lwd=2)
lines(xx1B[660:1080], yy1B[660:1080], lwd=2)
lines(xx3B[13:455], yy3B[13:455], lwd=2)
lines(xx3B[538:1055], yy3B[538:1055], lwd=2)
lines(xx3B[1079:1135], yy3B[1079:1135], lwd=2)

# Plot longitudes
vals <- seq(-7, 8, 1)*pi/8
for(lambda in vals){
  sl <- sin(lambda)
  cl <- cos(lambda)
  phi <- atan(((0-1)*k2*cl)/(k1))
  angs <- seq(phi, 5.6*pi/12, length=500)
  xx <- cos(angs)*sl
  yy <- k2*sin(angs)-k1*cos(angs)*cl
  lines(xx, yy, lwd=.5)
}


# Plot Latitudes
# vals2 <- seq(-2.8, 5.6, 1.4)*pi/12
vals2 <- c(-1.5, 0, 1.5, 3.0, 4.5, 5.6)*pi/12
for(phi in vals2){
  if (k1*sin(phi) > k2 * cos(phi)) 
    crit <- pi 
  else 
    crit <- acos((-k1*sin(phi))/(k2*cos(phi)))
  lambdas <- seq(-crit, crit, length=500)
  xx <- cos(phi)*sin(lambdas)
  yy <- k2*sin(phi)-k1 * cos(phi)*cos(lambdas)
  lines(xx, yy, lwd=.5)
}


# create lines for X1X2- and X1X3-planes
lambda <- pi/4
sl <- sin(lambda)
cl <- cos(lambda)
phi <- atan(((0-1)*k2*cl)/(k1))
angs <- seq(phi, pi+phi, length=500)
xx <- cos(angs)*sl
yy <- k2*sin(angs)-k1*cos(angs)*cl
lines(xx, yy, lwd=2)
lambda <- 3*pi/4
sl <- sin(lambda)
cl <- cos(lambda)
phi <- atan(((0-1)*k2*cl)/(k1))
angs <- seq(phi, pi+phi, length=500)
xx <- cos(angs)*sl
yy <- k2*sin(angs)-k1*cos(angs)*cl
lines(xx, yy, lwd=2)

# create line for X2X3-plane
phi <- 0
crit <- acos((-k1*sin(phi))/(k2 * cos(phi)))
lambdas <- seq(-crit, crit, length=500)
xx <- cos(phi)*sin(lambdas)
yy <- k2*sin(phi)-k1 * cos(phi)*cos(lambdas)
lines(xx, yy, lwd=2)

# create axes
lines(c(0.00, 0.00), c(k2*sin(pi/2), 1.11), lwd=4)
lines(c(0.00, 0.00), c(-1, -1.12), lwd=4)
a2x <- sin(-pi/4)
a2y <- cos(-pi/4)*(-k1)
lines(c(a2x, 1.5*a2x), c(a2y, 1.5*a2y), lwd=4)
a3x <- sin(pi/4)
a3y <- cos(pi/4)*(-k1)
lines(c(a3x, 1.5*a3x), c(a3y, 1.5*a3y), lwd=4)
k <- sqrt(a3x^2+a3y^2)
lines(c(-a3x/k, 1.2*(-a3x/k)), c(-a3y/k, 1.2*(-a3y/k)), lwd=4)
k <- sqrt(a2x^2+a2y^2)
lines(c(-a2x/k, 1.2*(-a2x/k)), c(-a2y/k, 1.2*(-a2y/k)), lwd=4)


# add text
text(-1.07, -.85, expression(bold(X[2])))
text(1.11, -.85, expression(bold(X[3])))
text(0.1, 1.12, expression(bold(X[1])))

lines(XX, YY, type="l")




}
figure1.12 <- function() {


#
# Comment:
#
# Code by Steve Miller 
# (Graduate student in the Statistics Department, The University of Auckland).
#
# An example of a one-off image using the traditional graphics system.  
# All parameters are hard-coded and the image only looks right when 
# drawn with a specific aspect ratio (4:1).
#
# Also an example of drawing an empty plot with a specific coordinate system
# and then building up a final image by drawing individual lines and
# and pieces of text.
# 
# Small point of interest is the use of some special glyphs (e.g., treble
# clef) from the Hershey vector fonts.
#



# TOP: music
par(yaxt = "n", xaxt = "n", ann = F, fig = c(0, 1, 0, 1), 
    mar = c(0, 0, 0, 0), cex=0.5)
plot(1:10, type = "n", xlab = "", ylab = "")
title(main = "A Little Culture", line = -1)
E = 5; F = 5.2; G = 5.4; A = 5.6; B = 5.8; C = 6; D = 6.2; E2 = 6.4; F2 = 6.6

# stave
for (i in c(E, G, B, D, F2)) {
	lines(x = c(1, 10), y = rep(i, 2))
}

# Hershey characters (treble clef, crotchet rest, sharp)
s1 = list(x = 1.2, y = G) #place clef on G
text(list(x = c(s1$x, s1$x + 8.5, s1$x + .5), y = c(s1$y, s1$y + .4, F2)), 
     vfont = c("serif", "plain"), 
     labels = c("\\#H2330", "\\#H2378", "\\#H2323"), 
     cex = 2) 

# time signature
text(x = rep(s1$x + .3, 2), y = c(s1$y, s1$y + .8), 
     labels = c("4", "4"), cex = 0.8)

# notes
points(list(y = c(B, A, G, A, B, B, B), 
            x = c(s1$x + 1, s1$x + 2, s1$x + 3, s1$x + 4, s1$x + 5.5, 
                  s1$x + 6.5, s1$x + 7.5)), 
       pch = 16, cex = 1.2)

# note tails
tail = 1.05
for (n in c(B, A, G, A)) {
  lines(x = rep(s1$x + tail, 2), y = c(n, n + 1))
  tail = tail + 1
}

tail = tail + .5
for (n in c(B, B, B)) {
  lines(x = rep(s1$x + tail, 2), y = c(n, n + 1))
  tail = tail + 1
}

# bar lines
lines(x = rep(1, 2), y = c(E, F2))
lines(x = rep(s1$x + 4.75, 2), y = c(E, F2))
lines(x = rep(9.9, 2), y = c(E, F2))
lines(x = rep(10, 2), y = c(E, F2), lwd = 2)

# lyrics
text(x = seq(s1$x + 1, s1$x + 8.5, by = 0.5), y = rep(4, 16), 
     labels = c("Ma-", "", "ry", "", "had", "", "a", "", "", 
                "lit-", "", "tle", "", "lamb", "", ""), 
     cex = 1, font = 4)



}
figure1.13 <- function() {


pic <- pixmap::read.pnm(system.file("extra", "AfterTheBombs.pnm", 
                                    package="RGraphics"))

source(system.file("extra", "as.raster.R", package="RGraphics"))

w <- 1024 # 578
h <- 768 # 500
picRaster <- as.raster(pic)
bg <- picRaster # [1:h, (1024 - w):1024]

unknown <- 8.7
total <- 9.1
known <- total - unknown

theta0 <- pi/4
thetaN <- theta0 + 2*pi*unknown/total
theta <- seq(theta0, thetaN, length.out=100)
x <- 0.3*c(0, cos(theta)) + 0.5
y <- 0.3*c(0, sin(theta)) + 0.35

# grid.newpage()
grid.raster(bg)
pushViewport(viewport(width=unit(1, "snpc"), height=unit(1, "snpc"),
                      gp=gpar(cex=1.2)))
grid.polygon(x, y, gp=gpar(col=NA, fill=rgb(.67, 0, .11, .7)))
label1 <- textGrob("UNACCOUNTED\nFOR",
                   unit(.2, "npc") - unit(2, "mm"),
                   unit(.6, "npc") + unit(2, "mm"),
                   gp=gpar(cex=1.4, fontface="bold"),
                   just=c("right", "bottom"))
grid.rect(.2, .6, just=c("right", "bottom"),
          width=grobWidth(label1) + unit(4, "mm"),
          height=grobHeight(label1) + unit(4, "mm"),
          gp=gpar(col=NA, fill=rgb(1, 1, 1, .5)))
grid.draw(label1)
label2 <- textGrob("ACCOUNTED\nFOR", 
                   unit(.8, "npc") + unit(2, "mm"),
                   unit(.6, "npc") + unit(2, "mm"),
                   gp=gpar(cex=1.4, fontface="bold"),
                   just=c("left", "bottom"))
grid.rect(.8, .6, just=c("left", "bottom"),
          width=grobWidth(label2) + unit(4, "mm"),
          height=grobHeight(label2) + unit(4, "mm"),
          gp=gpar(col=NA, fill=rgb(1, 1, 1, .5)))
grid.draw(label2)
grid.segments(c(.2, .8), .6,
              c(.3, .7), .5,
              gp=gpar(lty="dotted", lwd=2))
heading <- textGrob("The Department of Defense is unable to account for the use of
$8.7 billion of the $9.1 billion it spent on reconstruction in Iraq",
                    x=unit(0.5, "cm"),
                    y=unit(3, "lines"),
                    just=c("left", "top"),
                    gp=gpar(cex=1, col="white"))
pushViewport(viewport(x=0, y=1,
                      just=c("left", "top"),
                      height=grobHeight(heading) + unit(4, "lines"),
                      width=grobWidth(heading) + unit(1, "cm")))
grid.rect(gp=gpar(fill="black"))
grid.segments(x0=unit(0.5, "cm"),
              x1=unit(1, "npc") - unit(0.5, "cm"),
              y0=unit(1, "npc") - unit(2, "lines"),
              y1=unit(1, "npc") - unit(2, "lines"),
              gp=gpar(col="grey50", lwd=2))
grid.text("That's 96 Percent",
          x=unit(0.5, "cm"),
          y=unit(1, "npc") - unit(1, "lines"),
          just="left",
          gp=gpar(fontface="bold", col="white"))
grid.draw(heading)
popViewport(2)



}
figure11.1 <- function() {




par(mar=rep(1, 4))
par(mfrow=c(1, 2))
plot(faithful)
gplots::textplot(capture.output(summary(faithful)))




}
figure11.2 <- function() {





par(mar=rep(1, 4))
plot(pressure)
plotrix::addtable2plot(0, 300, pressure[13:19, ])




}
figure11.3 <- function() {




gridExtra::grid.table(pressure[13:19, ], show.box=TRUE, 
           separator="black")




}
figure11.4 <- function() {









par(mar=rep(1, 4))
x <- rnorm(20)
y <- rnorm(20)
plot(x, y, pch=16, col="gray")

xy <- plotrix::emptyspace(x, y)
text(xy, label="largest\nempty\nregion")

xy2 <- Hmisc::largest.empty(x, y, 1, 1)
rect(xy2$x - .5, xy2$y - .5, 
     xy2$x + .5, xy2$y + .5)




}
figure11.5 <- function() {
x <- runif(10)
y <- rnorm(10)




















par(mar=c(1, 1, 2, 1))
plot(x, y, pch=21, bg="gray", ylim=c(-3, 3), asp=1)
plotrix::spread.labels(x, y, labels=1:10)

mtext("spread.labels", side=3, line=0)



par(mar=c(1, 1, 2, 1))
plot(x, y, pch=21, bg="gray",  
     ylim=c(-2, 3), xlim=c(-.5, 1.5))
plotrix::thigmophobe.labels(x, y, labels=1:10)

mtext("thigmophobe.labels", side=3, line=0)



par(mar=c(1, 1, 2, 1))
plot(x, y, pch=21, bg="gray", ylim=c(-3, 3), asp=1)
adjy <- TeachingDemos::spread.labs(y, strheight("10", cex=1.5))
text(-0.5, adjy, labels=1:10, pos=2)
segments(-0.5, adjy, x, y)

mtext("spread.labs", side=3, line=0)



par(mar=c(1, 1, 2, 1))
plot(x, y, pch=16, col="gray", ylim=c(-2, 3), xlim=c(-.5, 1.5))
maptools::pointLabel(x, y, labels=as.character(1:10))

mtext("pointLabel", side=3, line=0)



par(mar=c(1, 1, 2, 1))
sx <- sort(x)
sy <- sort(y)
lines <- list(A=list(x=sx, y=y, lty=1), 
              B=list(x=sx, y=sy, lty=2),
              C=list(x=sx, y=rev(y), lty=3), 
              D=list(x=sx, y=rev(sy), lty=4))

plot(x, y, type="n", ylim=c(-3, 3))
lapply(lines, function(l) do.call("lines", l))
Hmisc::labcurve(lines)

mtext("labcurve", side=3, line=0)


}
figure11.6 <- function() {









par(mar=rep(0, 4))
plot.new()
plot.window(0:1, c(.1, 1))
plotrix::draw.circle(.1, .9, radius=1:5/100)
plotrix::draw.arc(.3, .9, radius=1:5/100, 
         deg1=45, deg2=seq(360, 160, -50))
plotrix::arctext("arctext", center=c(.5, .85), radius=.05,
        stretch=1.2)

text(.1, .8, "draw.circle")
text(.3, .8, "draw.arc")
plotrix::boxed.labels(.7, .85, "boxed.labels", bg="gray90")
plotrix::textbox(c(.85, 1), .9, "this is a textbox .")

plotrix::gradient.rect(.05, .5, .15, .7, col=gray(0:20/21))
plotrix::cylindrect(.25, .5, .35, .7, "black")
plotrix::rectFill(.45, .5, .55, .7, pch=16)

text(.1, .45, "gradient.rect")
text(.3, .45, "cylindrect")
text(.5, .45, "rectFill")
x <- c(.65, .65, .75, .75)
y <- c(.5, .7, .7, .5)
plotrix::polygon.shadow(x, y, offset=c(2/100, -2/100))
polygon(x, y, col="white")

text(.7, .45, "polygon.shadow")
TeachingDemos::shadowtext(.9, .6, "shadowtext")

TeachingDemos::my.symbols(seq(.3, .7, .2), .3,
           TeachingDemos::ms.male, inches=.2)
TeachingDemos::my.symbols(c(.4, .6), .3,
           TeachingDemos::ms.female, inches=.2)

text(.5, .2, "my.symbols")
box(col="gray")



}
figure11.7 <- function() {
# library(gridExtra)
# 
# 
# 
# gridExtra::grid.ellipse(x=1:6/7, y=rep(.8, 6), size=.1, 
#              default.units="npc", size.unit="npc", 
#              ar=1:6, angle=1:6*15/180*pi)
# grid.text("grid.ellipse", y=.7)
# gridExtra::grid.pattern(x=1:6/7, y=.5, width=unit(.1, "npc"),
#              height=unit(.1, "npc"), pattern=1:6,
#              motif.cex=.7, gp=gpar(fill="gray80"))
# 
# grid.text("grid.pattern", y=.4)
# gridExtra::grid.barbed(1:6/7, y=rep(c(.15, .25), 3), 
#             size=unit(.05, "snpc"), 
#             pch=21, gp=gpar(fill="gray"))
# 
# grid.text("grid.barbed", y=.1)
# grid.rect(gp=gpar(col="gray", fill=NA))
# 
# 
# 
}
figure11.8 <- function() {




gplots::plotmeans(mpg ~ cyl, mtcars, 
          barcol="black", n.label=FALSE, connect=FALSE)



}
figure11.9 <- function() {




grid.rect(1:10/11, .75, width=1/15, height=1/3,
          gp=gpar(col=NA,
            fill=colorspace::sequential_hcl(10, 0, 0, c(20, 90))))

grid.rect(1:10/11, .25, width=1/15, height=1/3,
          gp=gpar(col=NA,
            fill=colorspace::diverge_hcl(10, 0, 0, c(20, 90))))

grid.rect(gp=gpar(col="gray", fill=NA))



}
figure11.10 <- function() {









par(mar=rep(1, 4))
plot(rnorm(100), rnorm(100), pch=16, col="gray",
     ann=FALSE, axes=FALSE)
box()

plotrix::corner.label("top-left", x=-1, y=1)

gplots::smartlegend(x="right", y="top", 
            legend="top-right", pch=16, 
            col="gray", bg="white")

text(grconvertX(1, "npc"), grconvertY(0, "npc"), 
     adj=c(1, 0), labels="bottom-right")




}
figure11.11 <- function() {




plot(window(Nile, 1920, 1940))
TeachingDemos::subplot({ plot(Nile, axes=FALSE, ann=FALSE)
          rect(1920, 0, 1940, 2000, border=NA, col="gray")
          box()
          lines(Nile) }, 
        x=1920, y=1000, size=c(1.5, .75), hadj=0)



}
figure11.12 <- function() {
kelvin <- pressure$temperature + 273.15








with(pressure,
     {
         plot(temperature, pressure, axes=FALSE)
         axis(2)
         box()
         staxlab(1, at=temperature, cex=.7)
     })



with(pressure,
     {
         plot(kelvin, pressure, xlim=c(250, 650))
         axis.break(1)
     })


}
figure11.13 <- function() {
kelvin <- pressure$temperature + 273.15












with(pressure,
     revaxis(temperature, pressure))



plot(kelvin, pressure$pressure)
TeachingDemos::updateusr(c(0, 1), 0:1, c(-273.15, -272.15), 0:1)
abline(v=100)
text(x=100, y=700, " water boils", adj=0)



pdf("Figures/extra-zoomplot-%d.pdf", onefile=FALSE,
    width=4, height=4)
dev.control("enable")
plot(pressure)
TeachingDemos::zoomplot(c(0, 150), c(0, 3))

dev.off()
png("Web/extra-zoomplot%d.png", width=320, height=320)
dev.control("enable")
plot(pressure)
TeachingDemos::zoomplot(c(0, 150), c(0, 3))

dev.off()
system("cp Web/extra-zoomplot2.png Web/extra-axisscale3.png")



}
figure5.1 <- function() {
print(
qplot(temperature, pressure, data=pressure)

)


}
figure5.2 <- function() {
print(
qplot(temperature, pressure, data=pressure,
      main="Vapor Pressure of Mercury",
      geom=c("point", "line"))

)


}
figure5.3 <- function() {
# grid.newpage()
layvp <- viewport(layout=grid.layout(1, 5,
                    heights=unit(1, "inch"),
                    widths=unit(c(1, .5), "inch")),
                  name="vplay")
vpi <- function(i) {
    viewport(layout.pos.col=i, name=paste("vp", i, sep=""))
}
pushViewport(layvp)
pushViewport(viewport(layout.pos.col=1:5))
grid.rect(width=1.1, gp=gpar(col=NA, fill="gray"))
popViewport()
pushViewport(vpi(1))
upViewport()
pushViewport(vpi(3))
upViewport()
pushViewport(vpi(5))
upViewport(0)
for (i in c(1, 3, 5)) {
    grid.roundrect(height=.5,
                   vp=paste("vplay::vp", i, sep=""),
                   name=paste("rr", i, sep=""),
                   gp=gpar(fill="white"))
}
grid.text("data", vp="vplay::vp1")
grid.text("aesthetic", vp="vplay::vp3")
grid.text("geom", vp="vplay::vp5")
arr <- arrow(length=unit(3, "mm"), type="closed")
grid.segments(grobX("rr1", 0), .5,
              grobX("rr3", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr3", 0), .5,
              grobX("rr5", 180), .5,
              arrow=arr, gp=gpar(fill="black"))



}
figure5.4 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg))

)


print(
p + geom_point(aes(x=disp, y=mpg, shape=gear),
               size=4) +
    theme(legend.position="none")
)


print(
p + geom_text(aes(x=disp, y=mpg, label=gear))

)


lmcoef <- coef(lm(mpg ~ disp, mtcars2))



print(
p + geom_point(aes(x=disp, y=mpg)) +
    geom_abline(intercept=lmcoef[1], slope=lmcoef[2])

)


}
figure5.5 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg)) +
    scale_y_continuous(name="miles per gallon") +
    scale_x_continuous(name="displacement (cu.in.)")

)


print(
p + geom_point(aes(x=disp, y=mpg)) +
    scale_y_continuous(limits=c(0, 40)) 

)


print(
p + geom_point(aes(x=disp, y=mpg, 
                   color=trans), size=4) +
    scale_color_manual(values=c(automatic=gray(2/3),
                         manual=gray(1/3)))

)


}
figure5.6 <- function() {
# grid.newpage()
layvp <- viewport(layout=grid.layout(1, 7,
                    heights=unit(1, "inch"),
                    widths=unit(c(1, .5), "inch")),
                  name="vplay")
vpi <- function(i) {
    viewport(layout.pos.col=i, name=paste("vp", i, sep=""))
}
pushViewport(layvp)
pushViewport(viewport(layout.pos.col=1:7))
grid.rect(width=1.1, gp=gpar(col=NA, fill="gray"))
popViewport()
pushViewport(vpi(1))
upViewport()
pushViewport(vpi(3))
upViewport()
pushViewport(vpi(5))
upViewport()
pushViewport(vpi(7))
upViewport(0)
for (i in c(1, 3, 5, 7)) {
    grid.roundrect(height=.5,
                   vp=paste("vplay::vp", i, sep=""),
                   name=paste("rr", i, sep=""),
                   gp=gpar(fill="white"))
}
grid.text("data", vp="vplay::vp1")
grid.text("scale", vp="vplay::vp3")
grid.text("aesthetic", vp="vplay::vp5")
grid.text("geom", vp="vplay::vp7")
arr <- arrow(length=unit(3, "mm"), type="closed")
grid.segments(grobX("rr1", 0), .5,
              grobX("rr3", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr3", 0), .5,
              grobX("rr5", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
# grid.curve(grobX("rr1", 27), grobY("rr1", 27), 
#            grobX("rr5", 153), grobY("rr5", 153),
#            square=FALSE, ncp=8, curvature=-.3,
#            arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr5", 0), .5,
              grobX("rr7", 180), .5,
              arrow=arr, gp=gpar(fill="black"))



}
figure5.7 <- function() {
# grid.newpage()
layvp <- viewport(layout=grid.layout(1, 9,
                    heights=unit(1, "inch"),
                    widths=unit(c(1, .5), "inch")),
                  name="vplay")
vpi <- function(i) {
    viewport(layout.pos.col=i, name=paste("vp", i, sep=""))
}
pushViewport(layvp)
pushViewport(viewport(layout.pos.col=1:9))
grid.rect(width=1.1, gp=gpar(col=NA, fill="gray"))
popViewport()
pushViewport(vpi(1))
upViewport()
pushViewport(vpi(3))
upViewport()
pushViewport(vpi(5))
upViewport()
pushViewport(vpi(7))
upViewport()
pushViewport(vpi(9))
upViewport(0)
for (i in c(1, 3, 5, 7, 9)) {
    grid.roundrect(height=.5,
                   vp=paste("vplay::vp", i, sep=""),
                   name=paste("rr", i, sep=""),
                   gp=gpar(fill="white"))
}
grid.text("data", vp="vplay::vp1")
grid.text("scale", vp="vplay::vp3")
grid.text("stat", vp="vplay::vp5")
grid.text("aesthetic", vp="vplay::vp7")
grid.text("geom", vp="vplay::vp9")
arr <- arrow(length=unit(3, "mm"), type="closed")
grid.segments(grobX("rr1", 0), .5,
              grobX("rr3", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr3", 0), .5,
              grobX("rr5", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr5", 0), .5,
              grobX("rr7", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr7", 0), .5,
              grobX("rr9", 180), .5,
              arrow=arr, gp=gpar(fill="black"))



}
figure5.8 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_bar(aes(x=trans))

)


update_geom_defaults("smooth", aes(color="black"))
print(
p + geom_smooth(aes(x=disp, y=mpg))

)


}
figure5.9 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg, shape=trans)) +
    scale_shape_manual(values=c(1, 3))

)


print(
ggplot(mtcars2, aes(x=disp, y=mpg)) + 
    geom_point() +
    stat_smooth(aes(group=trans),
                method="lm")

)


}
figure5.10 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_bar(aes(x=trans, fill=factor(cyl)),
             color="black") +
    scale_fill_manual(values=gray(1:3/3))

)


print(
p + geom_bar(aes(x=trans, fill=factor(cyl)),
             color="black",
             position="dodge") +
    scale_fill_manual(values=gray(1:3/3))

)


print(
p + geom_bar(aes(x=trans, fill=factor(cyl)),
             color="black",
             position="fill") +
    scale_fill_manual(values=gray(1:3/3))

)


}
figure5.11 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg)) + 
    scale_x_continuous(trans="log") +
    scale_y_continuous(trans="log") +
    geom_line(aes(x=disp, y=mpg), stat="smooth", 
              method="lm")

)


print(
p + geom_point(aes(x=disp, y=mpg)) + 
    scale_x_continuous(trans="log") +
    scale_y_continuous(trans="log") +
    geom_line(aes(x=disp, y=mpg), stat="smooth", 
              method="lm") +
    coord_trans(xtrans="exp", ytrans="exp")

)


print(
p + geom_bar(aes(x="", fill=trans)) +
    scale_fill_manual(values=gray(1:2/3))

)


print(
p + geom_bar(aes(x="", fill=trans)) +
    scale_fill_manual(values=gray(1:2/3)) +
    coord_polar(theta="y") 

)


}
figure5.12 <- function() {
# grid.newpage()
layvp <- viewport(layout=grid.layout(1, 11,
                    heights=unit(1, "inch"),
                    widths=unit(c(1, .5), "inch")),
                  name="vplay")
vpi <- function(i) {
    viewport(layout.pos.col=i, name=paste("vp", i, sep=""))
}
pushViewport(layvp)
pushViewport(viewport(layout.pos.col=1:11))
grid.rect(width=1.1, gp=gpar(col=NA, fill="gray"))
popViewport()
pushViewport(vpi(1))
upViewport()
pushViewport(vpi(3))
upViewport()
pushViewport(vpi(5))
upViewport()
pushViewport(vpi(7))
upViewport()
pushViewport(vpi(9))
upViewport()
pushViewport(vpi(11))
upViewport(0)
for (i in c(1, 3, 5, 7, 9, 11)) {
    grid.roundrect(height=.5,
                   vp=paste("vplay::vp", i, sep=""),
                   name=paste("rr", i, sep=""),
                   gp=gpar(fill="white"))
}
grid.text("data", vp="vplay::vp1")
grid.text("scale", vp="vplay::vp3")
grid.text("stat", vp="vplay::vp5")
grid.text("aesthetic", vp="vplay::vp7")
grid.text("geom", vp="vplay::vp9")
grid.text("coord", vp="vplay::vp11")
arr <- arrow(length=unit(3, "mm"), type="closed")
grid.segments(grobX("rr1", 0), .5,
              grobX("rr3", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr3", 0), .5,
              grobX("rr5", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr5", 0), .5,
              grobX("rr7", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr7", 0), .5,
              grobX("rr9", 180), .5,
              arrow=arr, gp=gpar(fill="black"))
grid.segments(grobX("rr9", 0), .5,
              grobX("rr11", 180), .5,
              arrow=arr, gp=gpar(fill="black"))



}
figure5.13 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg)) +
    facet_wrap(~ gear, nrow=1)

)


}
figure5.14 <- function() {
# mtcars2 <- mtcars
# mtcars2$trans <- factor(mtcars$am, 
#                         levels=0:1, 
#                         labels=c("automatic", "manual"))
# mtcars2$gear <- as.factor(mtcars$gear)
# mtcars2$am <- NULL
# mtcars2$vs <- NULL
# mtcars2$drat <- NULL
# mtcars2$carb <- NULL
# mtcars2$wt <- NULL
# mtcars2$hp <- NULL
# mtcars2$qsec <- NULL
# 
# # To keep R CMD check happy
# mpg <- mtcars2$mpg
# 
# 
# 
# p <- ggplot(mtcars2)
# 
# 
# 
# print(
# p + geom_point(aes(x=disp, y=mpg)) +
#     theme_bw()
# 
# )
# 
# 
# print(
# p + geom_point(aes(x=disp, y=mpg)) +
#     opts(axis.title.y=theme_text(angle=0))
# 
# )
# 
# 
# print(
# p + geom_point(aes(x=disp, y=mpg)) +
#     opts(axis.title.y=theme_blank())
# 
# )
# 
# 
# print(
# p + geom_point(aes(x=disp, y=mpg)) +
#     opts(title="Vehicle Fuel Efficiency")
# 
# )
# 
# 
}
figure5.15 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$gear <- as.factor(mtcars$gear)
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL
mtcars2$wt <- NULL
mtcars2$hp <- NULL
mtcars2$qsec <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



p <- ggplot(mtcars2)



print(
p + geom_point(aes(x=disp, y=mpg)) +
    geom_hline(yintercept=29)

)


gcLimits <- 
    data.frame(category=c("2WD car",
                 "4WD car",
                 "2WD small pick-up truck",
                 "4WD small pick-up truck",
                 "2WD std pick-up truck",
                 "4WD std pick-up truck"),
               limit=c(29, 24, 20, 18, 17, 16))


print(
p + geom_point(aes(x=disp, y=mpg)) +
    geom_hline(data=gcLimits, 
               aes(yintercept=limit),
               linetype="dotted") +
    geom_text(data=gcLimits,
              aes(y=limit + .1, label=category),
              x=70, hjust=0, vjust=0, size=3)

)


}
figure15.1 <- function() {




nodes <- c("grDevices", "graphics", "grid",
           "lattice", "ggplot2")
edgeList <- 
    list(grDevices=list(edges=c("graphics", "grid")),
         graphics=list(),
         grid=list(edges=c("lattice", "ggplot2")),
         lattice=list(),
         ggplot2=list())
simpleGNEL <- new("graphNEL",
                  nodes=nodes,
                  edgeL=edgeList,
                  edgemode="directed")




# Weird stuff happening if don't pre-layout graph
temp <- Rgraphviz::agopen(simpleGNEL, "")
Rgraphviz::plot(temp)



}
# placeholder



figure15.3 <- function() {




nodes <- c("grDevices", "graphics", "grid",
           "lattice", "ggplot2")
edgeList <- 
    list(grDevices=list(edges=c("graphics", "grid")),
         graphics=list(),
         grid=list(edges=c("lattice", "ggplot2")),
         lattice=list(),
         ggplot2=list())
simpleGNEL <- new("graphNEL",
                  nodes=nodes,
                  edgeL=edgeList,
                  edgemode="directed")



# Weird stuff happening if don't pre-layout graph
tempGraph <- Rgraphviz::agopen(simpleGNEL, "", layoutType="neato")
Rgraphviz::plot(tempGraph)



}
figure15.4 <- function() {




nodes <- c("grDevices", "graphics", "grid",
           "lattice", "ggplot2")
edgeList <- 
    list(grDevices=list(edges=c("graphics", "grid")),
         graphics=list(),
         grid=list(edges=c("lattice", "ggplot2")),
         lattice=list(),
         ggplot2=list())
simpleGNEL <- new("graphNEL",
                  nodes=nodes,
                  edgeL=edgeList,
                  edgemode="directed")



Rgraphviz::plot(simpleGNEL, 
     edgeAttrs=list(lty=c(`grDevices~graphics`="solid", 
                      `grDevices~grid`="solid",
                      `grid~lattice`="dashed", 
                      `grid~ggplot2`="dashed")),
     nodeAttrs=list(fillcolor=c(grDevices="white", 
                      graphics="gray90", grid="gray90",
                      lattice="gray60", ggplot2="gray60")))



}
figure15.5 <- function() {




load(system.file("extra", "grd.rda", package="RGraphics"))
grDraw <- function(layout) {
    ragrd <- Rgraphviz::agopen(grd, "", layoutType=layout)
    xy <- Rgraphviz::getNodeXY(ragrd)
    grid.newpage()
    pushViewport(viewport(width=1.1, height=1.1),
                 plotViewport(xscale=range(xy$x), yscale=range(xy$y)))
    grid.circle(xy$x, xy$y, default.units="native", 
                r=unit(.25, "mm"), gp=gpar(fill="black"))
    grdNodes <- graph::nodes(grd)
    grdEdges <- graph::edges(grd)
    mapply(function(start, ends) {
        if (length(ends) > 0) {
            grid.segments(xy$x[grdNodes == start],
                          xy$y[grdNodes == start],
                          xy$x[grdNodes %in% ends],
                          xy$y[grdNodes %in% ends],
                          default.units="native",
                          gp=gpar(col=rgb(0,0,0,.5)))
        }
    },
           as.list(grdNodes),
           grdEdges)
    for (i in c("grDevices", "graphics", "grid", "lattice", "ggplot2")) {
        grid.rect(xy$x[grdNodes == i], xy$y[grdNodes == i],
                  width=stringWidth(i), height=unit(1, "lines"),
                  default.units="native",
                  gp=gpar(col=NA, fill=rgb(.5, .5, .5, .5)))
        grid.text(i, xy$x[grdNodes == i], xy$y[grdNodes == i],
                  default.units="native",
                  gp=gpar(col="white"))
    }
}

png("Figures/graph-pkgdep.png", width=1350, height=1350, res=300)
grDraw("neato")
dev.off()
system("cp Figures/graph-pkgdep.png Web/")



}
figure15.6 <- function() {




nodes <- c("grDevices", "graphics", "grid",
           "lattice", "ggplot2")
edgeList <- 
    list(grDevices=list(edges=c("graphics", "grid")),
         graphics=list(),
         grid=list(edges=c("lattice", "ggplot2")),
         lattice=list(),
         ggplot2=list())
simpleGNEL <- new("graphNEL",
                  nodes=nodes,
                  edgeL=edgeList,
                  edgemode="directed")



Rgraphviz::toFile(Rgraphviz::agopen(simpleGNEL, ""), 
       filename="Figures/graph-graphvizrender.ps", 
       fileType="ps")



}
figure15.7 <- function() {





dh <- hypergraph::DirectedHyperedge(c("A", "B"), c("C", "D"))
hg <- hypergraph::Hypergraph(LETTERS[1:4], list(dh))
getMethod("plot", "graphBPH")(hyperdraw::graphBPH(hg))



}
figure15.8 <- function() {




treeIgraph <- igraph::graph.tree(10)
fullIgraph <- igraph::graph.full(10)



# See ?igraph.plotting for useful graph attributes
treeIgraph <- igraph::set.vertex.attribute(treeIgraph, "color", value="black")
treeIgraph <- igraph::set.edge.attribute(treeIgraph, "color", value="black")
plot(treeIgraph, 
     layout=igraph::layout.reingold.tilford(treeIgraph, root=1, flip.y=FALSE))



fullIgraph <- igraph::set.vertex.attribute(fullIgraph, "color", value="black")
fullIgraph <- igraph::set.edge.attribute(fullIgraph, "color", value="black")
plot(fullIgraph, layout=igraph::layout.circle)



}
figure15.9 <- function() {








nodes <- c("grDevices", "graphics", "grid",
           "lattice", "ggplot2")
edgeList <- 
    list(grDevices=list(edges=c("graphics", "grid")),
         graphics=list(),
         grid=list(edges=c("lattice", "ggplot2")),
         lattice=list(),
         ggplot2=list())
simpleGNEL <- new("graphNEL",
                  nodes=nodes,
                  edgeL=edgeList,
                  edgemode="directed")



simpleNetwork <- 
    network::network(rbind(c(1, 2),
                  c(1, 3),
                  c(3, 4),
                  c(3, 5)),
            vertex.attr=list(vertex.names=nodes))




par(mar=rep(2, 4), xpd=NA)
set.seed(2500)
plot(simpleNetwork, mode="fruchtermanreingold", 
     vertex.col=1, displaylabels=TRUE)




}
figure15.10 <- function() {




par(mar=rep(1, 4))
plot.new()

nodePos <- diagram::coordinates(c(2, 2, 2, 2))

diagram::straightarrow(nodePos[1, ], nodePos[3,])

diagram::straightarrow(nodePos[3, ], nodePos[4,])
diagram::straightarrow(nodePos[3, ], nodePos[5,])
diagram::straightarrow(nodePos[5, ], nodePos[6,])
diagram::straightarrow(nodePos[6, ], nodePos[4,])
diagram::straightarrow(nodePos[5, ], nodePos[7,])
diagram::straightarrow(nodePos[6, ], nodePos[8,])
diagram::straightarrow(nodePos[7, ], nodePos[8,])

diagram::textplain(nodePos[3, ] + c(.2, .02), lab="yes")
diagram::textplain(nodePos[5, ] + c(.2, .02), lab="yes")
diagram::textplain(nodePos[6, ] + c(.03, .15), lab="yes")
diagram::textplain(nodePos[3, ] + c(.2, .02), lab="yes")
diagram::textplain(nodePos[5, ] + c(-.03, .125), lab="no")
diagram::textplain(nodePos[7, ] + c(-.03, .125), lab="no")
diagram::textplain(nodePos[7, ] + c(.2, -.02), lab="no")
diagram::textplain(nodePos[8, ] + c(-.03, .125), lab="no")

diagram::textrect(nodePos[1, ], .05, .025, lab="start")

diagram::textdiamond(nodePos[3, ], .15, .1)
diagram::textplain(nodePos[3, ], .08,
          lab=c("do you", "understand flow", "charts?"))
diagram::textellipse(nodePos[4, ], .08, .08,
            lab=c("let's go", "drink."))
diagram::textdiamond(nodePos[5, ], .15, .1)
diagram::textplain(nodePos[5, ], .08,
          lab=c("you see", "the lines labeled", "'yes'?"))
diagram::textdiamond(nodePos[6, ], .15, .1)
diagram::textplain(nodePos[6, ], .08,
          lab=c("you see", "the lines labeled", "'no'?"))
diagram::textdiamond(nodePos[7, ], .15, .1)
diagram::textplain(nodePos[7, ], .08,
          lab=c("you see", "the lines labeled", "'no'?"))
diagram::textellipse(nodePos[8, ], .07, .07,
            lab=c("I hate", "you."))



}
figure18.1 <- function() {



}
figure18.2 <- function() {




source(system.file("extra", "as.raster.R", package="RGraphics"))




moonPhase <- function(x, y, phase, size=.05) {
  # size is in inches
  n <- 17
  angle <- seq(0, 2*pi, length=n)
  xx <- x + cos(angle)*xinch(size)
  yy <- y + sin(angle)*yinch(size)
  if (phase == "New")
    fill <- "black"
  else
    fill <- "white"
  polygon(xx, yy, col=fill)
  if (phase == "1Q")
    polygon(xx[(n/4):(n*3/4) + 1],
            yy[(n/4):(n*3/4) + 1],
            col="black")
  if (phase == "3Q")
    polygon(xx[c(1:(n/4 + 1), (n*3/4 + 1):n)],
            yy[c(1:(n/4 + 1), (n*3/4 + 1):n)],
            col="black")
}


# Original image from NASA
# http://grin.hq.nasa.gov/ABSTRACTS/GPN-2000-000473.html
rasterMoon <- pixmap::read.pnm(system.file("extra", "GPN-2000-000473.pgm",
                                   package="RGraphics"))
par(pin=c(3.5, 1.75), oma=c(0, 3, 0, 0), xaxs="i", yaxs="i", cex=.7)
plot.new()
rect(0, 0, 1, 1, col="black")
rasterImage(rasterMoon, .25, 0, .75, 1*813/703)
par(new=TRUE, xaxs="r", yaxs="r", las=1)
plot(lowTideDate, lowTideHour, type="n",
     ylim=range(mainHours), axes=FALSE, ann=FALSE)
# dashed reference lines
abline(v=phases$date,
       col="white", lty="dashed")
for (subset in list(1:13, 14:29, 30:31)) {
  lines(lowTideDate[subset], lowTideHour[subset],
        lwd=2, col="white")
  points(lowTideDate[subset], lowTideHour[subset],
         pch=16, col="white")
}
box()
axis.POSIXct(1, lowTideDate)
axis.POSIXct(2, at=mainHours, format="%H:%M")
mtext("Time of Low Tide (NZDT)", side=2, line=4, las=0, cex=.7)
mtext("Auckland, New Zealand January 2010", side=1, line=3, cex=.7)
axis(3, at=phases$date, labels=FALSE)
par(xpd=NA)
ymax <- par("usr")[4]
for (i in 1:nrow(phases))
    moonPhase(phases$date[i], ymax + yinch(.2), 
              phases$phase[i])
mtext("Phases of the Moon", side=3, line=3, cex=.7)








grid.moonPhase <- function(x, y, phase, size=unit(.05, "in")) {
  n <- 17
  angle <- seq(0, 2*pi, length=n)
  xx <- x + cos(angle)*size
  yy <- y + sin(angle)*size
  if (phase == "New")
    fill <- "black"
  else
    fill <- "white"
  grid.polygon(xx, yy, gp=gpar(fill=fill))
  if (phase == "1Q")
      grid.polygon(xx[(n/4):(n*3/4) + 1],
                   yy[(n/4):(n*3/4) + 1],
                   gp=gpar(fill="black"))
  if (phase == "3Q")
      grid.polygon(xx[c(1:(n/4 + 1), (n*3/4 + 1):n)],
                   yy[c(1:(n/4 + 1), (n*3/4 + 1):n)],
                   gp=gpar(fill="black"))
}

# grid.newpage()
pushViewport(viewport(gp=gpar(cex=0.7)),
             plotViewport(c(4, 5, 3, 1)),
             dataViewport(as.numeric(lowTideDate), 
                          as.numeric(mainHours)))
vectorMoon <- 
    grImport::readPicture(system.file("extra", "comic_moon.ps.xml",
                            package="RGraphics"))
grImport::grid.picture(vectorMoon)
grid.segments(unit(phases$date, "native"), 0,
              unit(phases$date, "native"), 1,
              gp=gpar(lty="dashed"))
for (subset in list(1:13, 14:29, 30:31)) {
  grid.lines(lowTideDate[subset], lowTideHour[subset],
             default.units="native", 
             gp=gpar(lwd=2))
  grid.points(lowTideDate[subset], lowTideHour[subset],
              pch=16, size=unit(2, "mm"))
}
grid.rect(gp=gpar(fill=NA))
xTicks <- seq(min(lowTideDate), max(lowTideDate), by="week")
grid.xaxis(at=xTicks, label=format(xTicks, "%b %d"))
grid.yaxis(at=mainHours, label=format(mainHours, "%H:%M"))
grid.text("Time of Low Tide (NZDT)", 
          x=unit(-4, "lines"), rot=90)
grid.text("Auckland, New Zealand January 2010", 
          y=unit(-3, "lines"))
grid.xaxis(main=FALSE, at=phases$date, label=FALSE)
for (i in 1:nrow(phases))
    grid.moonPhase(unit(phases$date[i], "native"),
                   unit(1, "npc") + unit(1, "lines"), 
                   phases$phase[i])
grid.text("Phases of the Moon", 
          y=unit(1, "npc") + unit(2, "lines"))
popViewport(2)



}
figure18.3 <- function() {




source(system.file("extra", "as.raster.R", package="RGraphics"))




moon <- pixmap::read.pnm(system.file("extra", "GPN-2000-000473.pgm",
                                   package="RGraphics"))
helmet <- pixmap::read.pnm(system.file("extra", "astronaut.pgm",
                               package="RGraphics"))

moonMatrix <- as.matrix(as.raster(moon))
helmetMatrix <- as.matrix(as.raster(helmet))

moonCrop <- moonMatrix[120:(119 + nrow(helmetMatrix)),
                       10:(9 + ncol(helmetMatrix))]
moonGreys <- col2rgb(moonCrop)[1, ]
helmetRGB <- col2rgb(helmetMatrix)
helmetMask <- matrix(rgb(helmetRGB[1, ],
                         helmetRGB[2, ],
                         helmetRGB[3, ],
                         moonGreys, maxColorValue=255), ncol=ncol(helmetMatrix))




pushViewport(viewport(layout=grid.layout(1, 2, respect=TRUE)))
pushViewport(viewport(layout.pos.row=1, layout.pos.col=1),
             viewport(width=.8, height=.8))
grid.raster(helmet)
popViewport(2)
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2),
             viewport(width=.8, height=.8))
grid.raster(moonCrop)
popViewport(2)



}
figure18.4 <- function() {




source(system.file("extra", "as.raster.R", package="RGraphics"))




moon <- pixmap::read.pnm(system.file("extra", "GPN-2000-000473.pgm",
                                   package="RGraphics"))
helmet <- pixmap::read.pnm(system.file("extra", "astronaut.pgm",
                               package="RGraphics"))

moonMatrix <- as.matrix(as.raster(moon))
helmetMatrix <- as.matrix(as.raster(helmet))

moonCrop <- moonMatrix[120:(119 + nrow(helmetMatrix)),
                       10:(9 + ncol(helmetMatrix))]
moonGreys <- col2rgb(moonCrop)[1, ]
helmetRGB <- col2rgb(helmetMatrix)
helmetMask <- matrix(rgb(helmetRGB[1, ],
                         helmetRGB[2, ],
                         helmetRGB[3, ],
                         moonGreys, maxColorValue=255), ncol=ncol(helmetMatrix))




grid.rect(width=.99, height=.99)
grid.raster(helmetMask)



}
figure18.5 <- function() {





grImport::PostScriptTrace(system.file("extra", "comic_moon.ps",
                            package="RGraphics"))


vectorMoon <- grImport::readPicture("comic_moon.ps.xml")


grImport::picturePaths(vectorMoon[1:6], fill="white", 
             freeScales=TRUE, nr=2, nc=3)



}
figure18.6 <- function() {





grImport::PostScriptTrace(system.file("extra", "comic_moon.ps",
                            package="RGraphics"))


vectorMoon <- grImport::readPicture("comic_moon.ps.xml")


grImport::grid.picture(vectorMoon[1:4])



grImport::grid.picture(vectorMoon, use.gc=FALSE)



}
figure8.1 <- function() {




makeImageRect <- function(nrow, ncol, cols, byrow) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }  
  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name="image")
}

imageGrob <- function(nrow, ncol, cols, byrow=TRUE,
                       name=NULL, gp=NULL, vp=NULL) { 
  igt <- gTree(nrow=nrow, ncol=ncol, 
               cols=cols, byrow=byrow,
               children=gList(makeImageRect(nrow, ncol, 
                                            cols, byrow)),
               gp=gp, name=name, vp=vp, 
               cl="imageGrob") 
  igt
}

grid.imageGrob <- function(...) {
  igt <- imageGrob(...)
  grid.draw(igt)
}


makeOzViewports <- function(ozRegion) {
  vpStack(viewport(name="ozlay", layout=grid.layout(1, 1,
                     widths=diff(ozRegion$rangex),
                     heights=diff(ozRegion$rangey), 
                     respect=TRUE)),
          viewport(name="ozvp", layout.pos.row=1, 
                   layout.pos.col=1,
                   xscale=ozRegion$rangex, 
                   yscale=ozRegion$rangey, 
                   clip=TRUE))
}

makeOzLines <- function(ozRegion) {
  numLines <- length(ozRegion$lines)
  lines <- vector("list", numLines)
  index <- 1
  for(i in ozRegion$lines) {
    lines[[index]] <- linesGrob(i$x, i$y, 
                    default.units="native",
                    vp=vpPath("ozlay", "ozvp"), 
                    name=paste("ozlines", index, sep=""))
    index <- index + 1
  }
  do.call("gList", lines)
}

ozGrob <- function(ozRegion, name=NULL, gp=NULL, vp=NULL) {
  gTree(ozRegion=ozRegion, name=name, gp=gp, vp=vp, 
    childrenvp=makeOzViewports(ozRegion), 
    children=makeOzLines(ozRegion), 
    cl="ozGrob")
}

grid.ozGrob <- function(...) {
  grid.draw(ozGrob(...))
}


ozImage <- function(mapLong, mapLat, 
                    imageLong, imageLat, cols) {
  grob(mapLong=mapLong, mapLat=mapLat, 
       imageLong=imageLong, imageLat=imageLat, cols=cols,
       cl="ozImage")  
}

drawDetails.ozImage <- function(x, recording) { 
  grid.draw(ozGrob(oz::ozRegion(xlim=x$mapLong, 
                            ylim=x$mapLat))) 
  depth <- downViewport(vpPath("ozlay", "ozvp"))
  pushViewport(viewport(y=min(x$imageLat), 
                        height=diff(range(x$imageLat)), 
                        x=max(x$imageLong), 
                        width=diff(range(x$imageLong)),
                        default.units="native", 
                        just=c("right", "bottom")))
  grid.draw(imageGrob(50, 50, cols=x$col)) 
  popViewport()
  upViewport(depth)
} 


calcBreaks <- function(nlevels, breaks, scale) {
  if (is.null(breaks)) {
    seq(min(scale), max(scale), diff(scale)/nlevels)
  } else {
    breaks
  }
}

ribbonVps <- function(nlevels, breaks, margin, scale) {
  breaks <- format(signif(calcBreaks(nlevels, breaks, scale), 
                          3))
  vpTree(
    viewport(name="layout", layout=
      grid.layout(3, 4,
        widths=unit.c(margin, unit(1, "line"),
                      max(unit(0.8, "line") + 
                          stringWidth(breaks)), margin),
        heights=unit.c(margin, unit(1, "null"), margin))),
    vpList(viewport(layout.pos.col=2, layout.pos.row=2,
                    yscale=scale, name="ribbon"),
           viewport(layout.pos.col=3, layout.pos.row=2,
                    yscale=scale, name="labels")))
}

ribbonKids <- function(nlevels, breaks, cols, scale) {
  breaks <- calcBreaks(nlevels, breaks, scale)
  nb <- length(breaks)
  tickloc <- breaks[-c(1, nb)]
  gList(rectGrob(y=unit(breaks[-1], "native"), 
                 height=unit(diff(breaks), "native"),
                 just="top", gp=gpar(fill=cols),
                 vp=vpPath("layout", "ribbon")),
        segmentsGrob(x1=unit(0.5, "line"),
                     y0=unit(tickloc, "native"),
                     y1=unit(tickloc, "native"),
                     vp=vpPath("layout", "labels")),
        textGrob(x=unit(0.8, "line"),
                 y=unit(tickloc, "native"),
                 just="left", 
                 label=format(signif(tickloc, 3)),
                 vp=vpPath("layout", "labels")))
}


ribbonLegend <- function(nlevels=NULL, breaks=NULL, cols, 
                         scale=range(breaks), 
                         margin=unit(0.5, "line"), 
                         gp=NULL, vp=NULL, name=NULL) {
  gTree(
    nlevels=nlevels, breaks=breaks, cols=cols, scale=scale, 
    children=ribbonKids(nlevels, breaks, cols, scale),
    childrenvp=ribbonVps(nlevels, breaks, margin, scale),
    gp=gp, vp=vp, name=name, cl="ribbonLegend")
}

widthDetails.ribbonLegend <- function(x) { 
  sum(layout.widths(viewport.layout(x$childrenvp[[1]]))) 
} 


mapLong <- c(132, 136)
mapLat <- c(-35, -31.5)
imageLong <- range(RGraphics::fluoro.predict$x)
imageLat <- range(RGraphics::fluoro.predict$y)
zbreaks <- seq(min(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               max(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               length=10)
zcol <- cut(RGraphics::fluoro.predict$z, zbreaks,
            include.lowest=TRUE, labels=FALSE)
ozgrays <- gray(0.5 + 1:9/20)
imageCols <- ozgrays[zcol]



ozKey <- function(x, y, width, height, just, 
                  mapLong, mapLat) {
  gTree(childrenvp=viewport(name="ozkeyframe",
                            x=x, y=y, just=just,
                            width=width, height=height),
        children=gList(ozGrob(oz::ozRegion(), vp="ozkeyframe",
                              gp=gpar(lwd=0.1)),
                       rectGrob(x=mean(mapLong),
                                y=mean(mapLat),
                                width=abs(diff(mapLong)),
                                height=abs(diff(mapLat)),
                                default.units="native",
                                gp=gpar(lwd=1),
                                vp=vpPath("ozkeyframe",
                                          "ozlay", "ozvp"))))
}


ozimage <- ozImage(mapLong, mapLat, 
                   imageLong, imageLat, imageCols)



ribbonlegend <- ribbonLegend(breaks=zbreaks, 
                             cols=ozgrays, 
                             scale=range(zbreaks),
                             gp=gpar(cex=0.7))



ozkey <- ozKey(x=unit(1, "npc") - unit(1, "mm"),
               y=unit(1, "npc") - unit(1, "mm"),
               width=unit(3.5, "cm"),
               height=unit(2, "cm"),
               just=c("right", "top"),
               mapLong, mapLat)



grid.rect(gp=gpar(col="gray"))
fg <- frameGrob()
fg <- packGrob(fg, ozimage)
fg <- placeGrob(fg, ozkey)
fg <- packGrob(fg, ribbonlegend, "right")
grid.draw(fg)




}
grid.imageFun <- function(nrow, ncol, cols, 
                          byrow=TRUE) {
  x <- (1:ncol)/ncol
  y <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(x, nrow)
    top <- rep(y, each=ncol)
  } else {
    right <- rep(x, each=nrow)
    top <- rep(y, ncol)
  }
  grid.rect(x=right, y=top,  
    width=1/ncol, height=1/nrow, 
    just=c("right", "top"),
    gp=gpar(col=NA, fill=cols),
    name="image") 
}


figure8.3 <- function() {
grays <- gray(0.5 + (rep(1:4, 4) - rep(0:3, each=4))/10)



grid.imageFun <- function(nrow, ncol, cols, 
                          byrow=TRUE) {
  x <- (1:ncol)/ncol
  y <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(x, nrow)
    top <- rep(y, each=ncol)
  } else {
    right <- rep(x, each=nrow)
    top <- rep(y, ncol)
  }
  grid.rect(x=right, y=top,  
    width=1/ncol, height=1/nrow, 
    just=c("right", "top"),
    gp=gpar(col=NA, fill=cols),
    name="image") 
}


pushViewport(viewport(layout=grid.layout(3, 5, widths=c(1,8,2,8,1),
  heights=unit(c(1, 8, 1), c("null", "null", "line")))))
pushViewport(viewport(layout.pos.col=2, 
                      layout.pos.row=2))
grid.imageFun(4, 4, grays)

popViewport()
pushViewport(viewport(layout.pos.col=2, 
                      layout.pos.row=3))
grid.text("(a)", gp=gpar(cex=0.7))
popViewport()
pushViewport(viewport(layout.pos.col=4,
                      layout.pos.row=2))
grid.imageFun(4, 4, grays, byrow=FALSE)

popViewport()
pushViewport(viewport(layout.pos.col=4, 
                      layout.pos.row=3))
grid.text("(b)", gp=gpar(cex=0.7))
popViewport(2)



}
grid.ozFun <- function(ozRegion) {
  pushViewport( 
    viewport(name="ozlay", 
             layout=grid.layout(1,1,
                      widths=diff(ozRegion$rangex),
                      heights=diff(ozRegion$rangey), 
                      respect=TRUE)))
  pushViewport(viewport(name="ozvp", 
                        layout.pos.row=1, 
                        layout.pos.col=1,
                        xscale=ozRegion$rangex, 
                        yscale=ozRegion$rangey, 
                        clip=TRUE)) 
  index <- 1
  for(i in ozRegion$lines) {
    grid.lines(i$x, i$y, default.units="native",
               name=paste("ozlines", index, sep="")) 
    index <- index + 1
  }
  upViewport(2) 
}


figure8.5 <- function() {




grid.ozFun <- function(ozRegion) {
  pushViewport( 
    viewport(name="ozlay", 
             layout=grid.layout(1,1,
                      widths=diff(ozRegion$rangex),
                      heights=diff(ozRegion$rangey), 
                      respect=TRUE)))
  pushViewport(viewport(name="ozvp", 
                        layout.pos.row=1, 
                        layout.pos.col=1,
                        xscale=ozRegion$rangex, 
                        yscale=ozRegion$rangey, 
                        clip=TRUE)) 
  index <- 1
  for(i in ozRegion$lines) {
    grid.lines(i$x, i$y, default.units="native",
               name=paste("ozlines", index, sep="")) 
    index <- index + 1
  }
  upViewport(2) 
}


grid.rect(gp=gpar(col="gray"))
grid.ozFun(oz::ozRegion())




}
figure8.6 <- function() {




grid.imageFun <- function(nrow, ncol, cols, 
                          byrow=TRUE) {
  x <- (1:ncol)/ncol
  y <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(x, nrow)
    top <- rep(y, each=ncol)
  } else {
    right <- rep(x, each=nrow)
    top <- rep(y, ncol)
  }
  grid.rect(x=right, y=top,  
    width=1/ncol, height=1/nrow, 
    just=c("right", "top"),
    gp=gpar(col=NA, fill=cols),
    name="image") 
}


grid.ozFun <- function(ozRegion) {
  pushViewport( 
    viewport(name="ozlay", 
             layout=grid.layout(1,1,
                      widths=diff(ozRegion$rangex),
                      heights=diff(ozRegion$rangey), 
                      respect=TRUE)))
  pushViewport(viewport(name="ozvp", 
                        layout.pos.row=1, 
                        layout.pos.col=1,
                        xscale=ozRegion$rangex, 
                        yscale=ozRegion$rangey, 
                        clip=TRUE)) 
  index <- 1
  for(i in ozRegion$lines) {
    grid.lines(i$x, i$y, default.units="native",
               name=paste("ozlines", index, sep="")) 
    index <- index + 1
  }
  upViewport(2) 
}


mapLong <- c(132, 136)
mapLat <- c(-35, -31.5)
imageLong <- range(RGraphics::fluoro.predict$x)
imageLat <- range(RGraphics::fluoro.predict$y)
zbreaks <- seq(min(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               max(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               length=10)
zcol <- cut(RGraphics::fluoro.predict$z, zbreaks,
            include.lowest=TRUE, labels=FALSE)
ozgrays <- gray(0.5 + 1:9/20)
imageCols <- ozgrays[zcol]



mapLong <- c(132, 136)
mapLat <- c(-35, -31.5)
imageLong <- range(RGraphics::fluoro.predict$x)
imageLat <- range(RGraphics::fluoro.predict$y)
zbreaks <- seq(min(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               max(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               length=10)
zcol <- cut(RGraphics::fluoro.predict$z, zbreaks,
            include.lowest=TRUE, labels=FALSE)
ozgrays <- gray(0.5 + 1:9/20)
imageCols <- ozgrays[zcol]

grid.rect(gp=gpar(col="gray"))
grid.ozFun(oz::ozRegion(xlim=mapLong, ylim=mapLat))

downViewport("ozvp")

pushViewport(viewport(y=min(imageLat), 
                      height=abs(diff(imageLat)), 
                      x=max(imageLong), 
                      width=abs(diff(imageLong)),
                      default.units="native", 
                      just=c("right", "bottom")))
grid.imageFun(50, 50, col=imageCols)
upViewport(0)




}
figure8.7 <- function() {




grid.imageFun <- function(nrow, ncol, cols, 
                          byrow=TRUE) {
  x <- (1:ncol)/ncol
  y <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(x, nrow)
    top <- rep(y, each=ncol)
  } else {
    right <- rep(x, each=nrow)
    top <- rep(y, ncol)
  }
  grid.rect(x=right, y=top,  
    width=1/ncol, height=1/nrow, 
    just=c("right", "top"),
    gp=gpar(col=NA, fill=cols),
    name="image") 
}


grid.ozFun <- function(ozRegion) {
  pushViewport( 
    viewport(name="ozlay", 
             layout=grid.layout(1,1,
                      widths=diff(ozRegion$rangex),
                      heights=diff(ozRegion$rangey), 
                      respect=TRUE)))
  pushViewport(viewport(name="ozvp", 
                        layout.pos.row=1, 
                        layout.pos.col=1,
                        xscale=ozRegion$rangex, 
                        yscale=ozRegion$rangey, 
                        clip=TRUE)) 
  index <- 1
  for(i in ozRegion$lines) {
    grid.lines(i$x, i$y, default.units="native",
               name=paste("ozlines", index, sep="")) 
    index <- index + 1
  }
  upViewport(2) 
}


mapLong <- c(132, 136)
mapLat <- c(-35, -31.5)
imageLong <- range(RGraphics::fluoro.predict$x)
imageLat <- range(RGraphics::fluoro.predict$y)
zbreaks <- seq(min(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               max(RGraphics::fluoro.predict$z, na.rm=TRUE), 
               length=10)
zcol <- cut(RGraphics::fluoro.predict$z, zbreaks,
            include.lowest=TRUE, labels=FALSE)
ozgrays <- gray(0.5 + 1:9/20)
imageCols <- ozgrays[zcol]



grid.rect(gp=gpar(col="gray"))
grid.ozFun(oz::ozRegion(xlim=mapLong, ylim=mapLat))

downViewport("ozvp")

pushViewport(viewport(y=min(imageLat), 
                      height=abs(diff(imageLat)), 
                      x=max(imageLong), 
                      width=abs(diff(imageLong)),
                      default.units="native", 
                      just=c("right", "bottom")))
grid.imageFun(50, 50, col=imageCols)
upViewport(0)

grid.edit("image", gp=gpar(fill=rev(ozgrays)[zcol]))
grid.gedit("^ozlines[0-9]+$", gp=gpar(col="gray", lwd=2))




}
makeImageRect <- function(nrow, ncol, cols, byrow) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }  
  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name="image")
}

imageGrob <- function(nrow, ncol, cols, byrow=TRUE,
                       name=NULL, gp=NULL, vp=NULL) { 
  igt <- gTree(nrow=nrow, ncol=ncol, 
               cols=cols, byrow=byrow,
               children=gList(makeImageRect(nrow, ncol, 
                                            cols, byrow)),
               gp=gp, name=name, vp=vp, 
               cl="imageGrob") 
  igt
}

grid.imageGrob <- function(...) {
  igt <- imageGrob(...)
  grid.draw(igt)
}


validDetails.imageGrob <- function(x) { 
  if (!is.numeric(x$nrow) || length(x$nrow) > 1 || 
      !is.numeric(x$ncol) || length(x$ncol) > 1)
    stop("nrow and ncol must be numeric and length 1")
  if (!is.logical(x$byrow))
    stop("byrow must be logical")
  x 
} 

validDetails.ozGrob <- function(x) {
  if (!inherits(x$ozRegion, "ozRegion"))
    stop("Invalid ozRegion")
  x
}


makeOzViewports <- function(ozRegion) {
  vpStack(viewport(name="ozlay", layout=grid.layout(1, 1,
                     widths=diff(ozRegion$rangex),
                     heights=diff(ozRegion$rangey), 
                     respect=TRUE)),
          viewport(name="ozvp", layout.pos.row=1, 
                   layout.pos.col=1,
                   xscale=ozRegion$rangex, 
                   yscale=ozRegion$rangey, 
                   clip=TRUE))
}

makeOzLines <- function(ozRegion) {
  numLines <- length(ozRegion$lines)
  lines <- vector("list", numLines)
  index <- 1
  for(i in ozRegion$lines) {
    lines[[index]] <- linesGrob(i$x, i$y, 
                    default.units="native",
                    vp=vpPath("ozlay", "ozvp"), 
                    name=paste("ozlines", index, sep=""))
    index <- index + 1
  }
  do.call("gList", lines)
}

ozGrob <- function(ozRegion, name=NULL, gp=NULL, vp=NULL) {
  gTree(ozRegion=ozRegion, name=name, gp=gp, vp=vp, 
    childrenvp=makeOzViewports(ozRegion), 
    children=makeOzLines(ozRegion), 
    cl="ozGrob")
}

grid.ozGrob <- function(...) {
  grid.draw(ozGrob(...))
}


ozImage <- function(mapLong, mapLat, 
                    imageLong, imageLat, cols) {
  grob(mapLong=mapLong, mapLat=mapLat, 
       imageLong=imageLong, imageLat=imageLat, cols=cols,
       cl="ozImage")  
}

drawDetails.ozImage <- function(x, recording) { 
  grid.draw(ozGrob(oz::ozRegion(xlim=x$mapLong, 
                            ylim=x$mapLat))) 
  depth <- downViewport(vpPath("ozlay", "ozvp"))
  pushViewport(viewport(y=min(x$imageLat), 
                        height=diff(range(x$imageLat)), 
                        x=max(x$imageLong), 
                        width=diff(range(x$imageLong)),
                        default.units="native", 
                        just=c("right", "bottom")))
  grid.draw(imageGrob(50, 50, cols=x$col)) 
  popViewport()
  upViewport(depth)
} 


editDetails.imageGrob <- function(x, specs) { 
  if (any(c("ncol", "nrow", "byrow") %in% names(specs))) { 
    x <- addGrob(x, makeImageRect(x$nrow, x$ncol,
                                  x$cols, x$byrow))
  } 
  if (any(c("cols") %in% names(specs))) { 
    x <- editGrob(x, "image", gp=gpar(fill=x$cols))
  } 
  x 
} 

editDetails.ozGrob <- function(x, specs) {
  if ("ozRegion" %in% names(specs)) {
    x$childrenvp <- makeOzViewports(x$ozRegion)
    x <- setChildren(x, makeOzLines(x$ozRegion))
  }
  x
}


figure8.13 <- function() {
grays <- gray(0.5 + (rep(1:4, 4) - rep(0:3, each=4))/10)



makeImageRect <- function(nrow, ncol, cols, byrow) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }  
  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name="image")
}

imageGrob <- function(nrow, ncol, cols, byrow=TRUE,
                       name=NULL, gp=NULL, vp=NULL) { 
  igt <- gTree(nrow=nrow, ncol=ncol, 
               cols=cols, byrow=byrow,
               children=gList(makeImageRect(nrow, ncol, 
                                            cols, byrow)),
               gp=gp, name=name, vp=vp, 
               cl="imageGrob") 
  igt
}

grid.imageGrob <- function(...) {
  igt <- imageGrob(...)
  grid.draw(igt)
}


editDetails.imageGrob <- function(x, specs) { 
  if (any(c("ncol", "nrow", "byrow") %in% names(specs))) { 
    x <- addGrob(x, makeImageRect(x$nrow, x$ncol,
                                  x$cols, x$byrow))
  } 
  if (any(c("cols") %in% names(specs))) { 
    x <- editGrob(x, "image", gp=gpar(fill=x$cols))
  } 
  x 
} 

editDetails.ozGrob <- function(x, specs) {
  if ("ozRegion" %in% names(specs)) {
    x$childrenvp <- makeOzViewports(x$ozRegion)
    x <- setChildren(x, makeOzLines(x$ozRegion))
  }
  x
}


pushViewport(viewport(layout=grid.layout(2, 1, 
                                         heights=unit(c(1, 1),
                                                      c("null", "line")),
                                         respect=TRUE)))
pushViewport(viewport(layout.pos.row=1))
grid.imageGrob(4, 4, grays, name="imageGrob")

popViewport()
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(a)", name="label")
popViewport()
grid.edit("imageGrob", byrow=FALSE)

grid.remove("label")
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(b)", name="label")
popViewport()
grid.edit("imageGrob::image", gp=gpar(col="white", lwd=6))

grid.remove("label")
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(c)", name="label")
popViewport(2)



}
figure8.14 <- function() {
grays <- gray(0.5 + (rep(1:4, 4) - rep(0:3, each=4))/10)



makeImageRect <- function(nrow, ncol, cols, byrow) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }  
  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name="image")
}

imageGrob <- function(nrow, ncol, cols, byrow=TRUE,
                       name=NULL, gp=NULL, vp=NULL) { 
  igt <- gTree(nrow=nrow, ncol=ncol, 
               cols=cols, byrow=byrow,
               children=gList(makeImageRect(nrow, ncol, 
                                            cols, byrow)),
               gp=gp, name=name, vp=vp, 
               cl="imageGrob") 
  igt
}

grid.imageGrob <- function(...) {
  igt <- imageGrob(...)
  grid.draw(igt)
}


editDetails.imageGrob <- function(x, specs) { 
  if (any(c("ncol", "nrow", "byrow") %in% names(specs))) { 
    x <- addGrob(x, makeImageRect(x$nrow, x$ncol,
                                  x$cols, x$byrow))
  } 
  if (any(c("cols") %in% names(specs))) { 
    x <- editGrob(x, "image", gp=gpar(fill=x$cols))
  } 
  x 
} 

editDetails.ozGrob <- function(x, specs) {
  if ("ozRegion" %in% names(specs)) {
    x$childrenvp <- makeOzViewports(x$ozRegion)
    x <- setChildren(x, makeOzLines(x$ozRegion))
  }
  x
}


pushViewport(viewport(layout=grid.layout(2, 1, 
                                         heights=unit(c(1, 1),
                                                      c("null", "line")),
                                         respect=TRUE)))
pushViewport(viewport(layout.pos.row=1))
grid.imageGrob(4, 4, grays, name="imageGrob")

grid.edit("imageGrob::image", gp=gpar(col="white"))

popViewport()
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(a)", name="label")
popViewport()
grid.edit("imageGrob", cols=rev(grays))

grid.remove("label")
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(b)", name="label")
popViewport()
grid.edit("imageGrob", byrow=FALSE)

grid.remove("label")
pushViewport(viewport(layout.pos.row=2, gp=gpar(cex=0.7)))
grid.text("(c)", name="label")
popViewport(2)



}
calcBreaks <- function(nlevels, breaks, scale) {
  if (is.null(breaks)) {
    seq(min(scale), max(scale), diff(scale)/nlevels)
  } else {
    breaks
  }
}

ribbonVps <- function(nlevels, breaks, margin, scale) {
  breaks <- format(signif(calcBreaks(nlevels, breaks, scale), 
                          3))
  vpTree(
    viewport(name="layout", layout=
      grid.layout(3, 4,
        widths=unit.c(margin, unit(1, "line"),
                      max(unit(0.8, "line") + 
                          stringWidth(breaks)), margin),
        heights=unit.c(margin, unit(1, "null"), margin))),
    vpList(viewport(layout.pos.col=2, layout.pos.row=2,
                    yscale=scale, name="ribbon"),
           viewport(layout.pos.col=3, layout.pos.row=2,
                    yscale=scale, name="labels")))
}

ribbonKids <- function(nlevels, breaks, cols, scale) {
  breaks <- calcBreaks(nlevels, breaks, scale)
  nb <- length(breaks)
  tickloc <- breaks[-c(1, nb)]
  gList(rectGrob(y=unit(breaks[-1], "native"), 
                 height=unit(diff(breaks), "native"),
                 just="top", gp=gpar(fill=cols),
                 vp=vpPath("layout", "ribbon")),
        segmentsGrob(x1=unit(0.5, "line"),
                     y0=unit(tickloc, "native"),
                     y1=unit(tickloc, "native"),
                     vp=vpPath("layout", "labels")),
        textGrob(x=unit(0.8, "line"),
                 y=unit(tickloc, "native"),
                 just="left", 
                 label=format(signif(tickloc, 3)),
                 vp=vpPath("layout", "labels")))
}


ribbonLegend <- function(nlevels=NULL, breaks=NULL, cols, 
                         scale=range(breaks), 
                         margin=unit(0.5, "line"), 
                         gp=NULL, vp=NULL, name=NULL) {
  gTree(
    nlevels=nlevels, breaks=breaks, cols=cols, scale=scale, 
    children=ribbonKids(nlevels, breaks, cols, scale),
    childrenvp=ribbonVps(nlevels, breaks, margin, scale),
    gp=gp, vp=vp, name=name, cl="ribbonLegend")
}

widthDetails.ribbonLegend <- function(x) { 
  sum(layout.widths(viewport.layout(x$childrenvp[[1]]))) 
} 


ozKey <- function(x, y, width, height, just, 
                  mapLong, mapLat) {
  gTree(childrenvp=viewport(name="ozkeyframe",
                            x=x, y=y, just=just,
                            width=width, height=height),
        children=gList(ozGrob(oz::ozRegion(), vp="ozkeyframe",
                              gp=gpar(lwd=0.1)),
                       rectGrob(x=mean(mapLong),
                                y=mean(mapLat),
                                width=abs(diff(mapLong)),
                                height=abs(diff(mapLat)),
                                default.units="native",
                                gp=gpar(lwd=1),
                                vp=vpPath("ozkeyframe",
                                          "ozlay", "ozvp"))))
}


figure8.18 <- function() {




makeOzViewports <- function(ozRegion) {
  vpStack(viewport(name="ozlay", layout=grid.layout(1, 1,
                     widths=diff(ozRegion$rangex),
                     heights=diff(ozRegion$rangey), 
                     respect=TRUE)),
          viewport(name="ozvp", layout.pos.row=1, 
                   layout.pos.col=1,
                   xscale=ozRegion$rangex, 
                   yscale=ozRegion$rangey, 
                   clip=TRUE))
}

makeOzLines <- function(ozRegion) {
  numLines <- length(ozRegion$lines)
  lines <- vector("list", numLines)
  index <- 1
  for(i in ozRegion$lines) {
    lines[[index]] <- linesGrob(i$x, i$y, 
                    default.units="native",
                    vp=vpPath("ozlay", "ozvp"), 
                    name=paste("ozlines", index, sep=""))
    index <- index + 1
  }
  do.call("gList", lines)
}

ozGrob <- function(ozRegion, name=NULL, gp=NULL, vp=NULL) {
  gTree(ozRegion=ozRegion, name=name, gp=gp, vp=vp, 
    childrenvp=makeOzViewports(ozRegion), 
    children=makeOzLines(ozRegion), 
    cl="ozGrob")
}

grid.ozGrob <- function(...) {
  grid.draw(ozGrob(...))
}


calcBreaks <- function(nlevels, breaks, scale) {
  if (is.null(breaks)) {
    seq(min(scale), max(scale), diff(scale)/nlevels)
  } else {
    breaks
  }
}

ribbonVps <- function(nlevels, breaks, margin, scale) {
  breaks <- format(signif(calcBreaks(nlevels, breaks, scale), 
                          3))
  vpTree(
    viewport(name="layout", layout=
      grid.layout(3, 4,
        widths=unit.c(margin, unit(1, "line"),
                      max(unit(0.8, "line") + 
                          stringWidth(breaks)), margin),
        heights=unit.c(margin, unit(1, "null"), margin))),
    vpList(viewport(layout.pos.col=2, layout.pos.row=2,
                    yscale=scale, name="ribbon"),
           viewport(layout.pos.col=3, layout.pos.row=2,
                    yscale=scale, name="labels")))
}

ribbonKids <- function(nlevels, breaks, cols, scale) {
  breaks <- calcBreaks(nlevels, breaks, scale)
  nb <- length(breaks)
  tickloc <- breaks[-c(1, nb)]
  gList(rectGrob(y=unit(breaks[-1], "native"), 
                 height=unit(diff(breaks), "native"),
                 just="top", gp=gpar(fill=cols),
                 vp=vpPath("layout", "ribbon")),
        segmentsGrob(x1=unit(0.5, "line"),
                     y0=unit(tickloc, "native"),
                     y1=unit(tickloc, "native"),
                     vp=vpPath("layout", "labels")),
        textGrob(x=unit(0.8, "line"),
                 y=unit(tickloc, "native"),
                 just="left", 
                 label=format(signif(tickloc, 3)),
                 vp=vpPath("layout", "labels")))
}


ribbonLegend <- function(nlevels=NULL, breaks=NULL, cols, 
                         scale=range(breaks), 
                         margin=unit(0.5, "line"), 
                         gp=NULL, vp=NULL, name=NULL) {
  gTree(
    nlevels=nlevels, breaks=breaks, cols=cols, scale=scale, 
    children=ribbonKids(nlevels, breaks, cols, scale),
    childrenvp=ribbonVps(nlevels, breaks, margin, scale),
    gp=gp, vp=vp, name=name, cl="ribbonLegend")
}

widthDetails.ribbonLegend <- function(x) { 
  sum(layout.widths(viewport.layout(x$childrenvp[[1]]))) 
} 


grid.ozGrob(oz::ozRegion())
downViewport("ozvp")
for (i in 1:(dim(RGraphics::ozTemp)[1])) {
  grid.points(RGraphics::ozTemp$long[i], RGraphics::ozTemp$lat[i], pch=16)
  rl <- ribbonLegend(breaks=c(min(RGraphics::ozTemp$min), 
                              RGraphics::ozTemp$min[i], 
                              RGraphics::ozTemp$max[i], 
                              max(RGraphics::ozTemp$max)),
                     cols=c("white", "gray", "white"),
                     gp=gpar(cex=.7))
  pushViewport(viewport(x=unit(RGraphics::ozTemp$long[i], "native"),
                        y=unit(RGraphics::ozTemp$lat[i], "native"),
                        height=unit(1, "in"),
                        width=grobWidth(rl),
                        clip="off"))
  grid.circle(r=1, 
              gp=gpar(col="gray", fill="white", alpha=0.8))
  grid.draw(rl)
  popViewport()
}
upViewport(0)



}
splitString <- function(text) {
  strings <- strsplit(text, " ")[[1]]
  newstring <- strings[1]
  linewidth <- stringWidth(newstring)
  gapwidth <- stringWidth(" ")
  availwidth <- 
    convertWidth(unit(1, "npc"), 
                 "in", valueOnly=TRUE) 
  for (i in 2:length(strings)) {
    width <- stringWidth(strings[i])
    if (convertWidth(linewidth + gapwidth + width, 
                     "in", valueOnly=TRUE) <
        availwidth) {
      sep <- " "
      linewidth <- linewidth + gapwidth + width
    } else {
      sep <- "\n"
      linewidth <- width
    }
    newstring <- paste(newstring, strings[i], sep=sep)
  }
  newstring
}   


figure8.20 <- function() {
splitString <- function(text) {
  strings <- strsplit(text, " ")[[1]]
  newstring <- strings[1]
  linewidth <- stringWidth(newstring)
  gapwidth <- stringWidth(" ")
  availwidth <- 
    convertWidth(unit(1, "npc"), 
                 "in", valueOnly=TRUE) 
  for (i in 2:length(strings)) {
    width <- stringWidth(strings[i])
    if (convertWidth(linewidth + gapwidth + width, 
                     "in", valueOnly=TRUE) <
        availwidth) {
      sep <- " "
      linewidth <- linewidth + gapwidth + width
    } else {
      sep <- "\n"
      linewidth <- width
    }
    newstring <- paste(newstring, strings[i], sep=sep)
  }
  newstring
}   


text <- "The quick brown fox jumps over the lazy dog."
grid.text(splitString(text), 
          x=0, y=1, just=c("left", "top")) 



splitTextGrob <- function(text, ...) {
  grob(text=text, cl="splitText", ...)
}

drawDetails.splitText <- function(x, recording) {
  grid.text(splitString(x$text),
            x=0, y=1, just=c("left", "top")) 
}


pushViewport(viewport(layout=grid.layout(2, 2)))
pushViewport(viewport(layout.pos.col=1))
pushViewport(viewport(width=0.5, height=0.9))
grid.rect(gp=gpar(col="gray"))
text <- "The quick brown fox jumps over the lazy dog."
grid.text(splitString(text), 
          x=0, y=1, just=c("left", "top")) 

popViewport(2)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
pushViewport(viewport(height=0.8))
grid.rect(gp=gpar(col="gray"))
splitText <- splitTextGrob(text, name="splitText")
grid.draw(splitText)

popViewport(2)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
pushViewport(viewport(height=0.8))
grid.rect(gp=gpar(col="gray"))
grid.draw(editGrob(splitText, gp=gpar(cex=1.5)))
popViewport(2)
popViewport()



}
splitTextGrob <- function(text, ...) {
  grob(text=text, cl="splitText", ...)
}

drawDetails.splitText <- function(x, recording) {
  grid.text(splitString(x$text),
            x=0, y=1, just=c("left", "top")) 
}


figure8.22 <- function() {
faceA <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.rect()
  grid.circle(x=c(0.25, 0.75), y=0.75, r=0.1)
  grid.lines(x=c(0.33, 0.67), y=0.25)
  popViewport()
}

faceB <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.draw(rectGrob())
  grid.draw(circleGrob(x=c(0.25, 0.75), y=0.75, r=0.1))
  grid.draw(linesGrob(x=c(0.33, 0.67), y=0.25))
  popViewport()
}


faceA(.5, .5, width=.1, height=.1)
angle <- seq(0, 2*pi, length=9)[-9]
for (i in angle) {
  x <- 0.5 + 0.3*cos(i)
  y <- 0.5 + 0.3*sin(i)
  faceA(x, y, 0.2*x, 0.2*y)
}
grid.rect(width=.9, height=.9, gp=gpar(col="gray", fill=NA))



}
faceA <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.rect()
  grid.circle(x=c(0.25, 0.75), y=0.75, r=0.1)
  grid.lines(x=c(0.33, 0.67), y=0.25)
  popViewport()
}

faceB <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.draw(rectGrob())
  grid.draw(circleGrob(x=c(0.25, 0.75), y=0.75, r=0.1))
  grid.draw(linesGrob(x=c(0.33, 0.67), y=0.25))
  popViewport()
}


faceC <- function(x, y, width, height) {
  gTree(childrenvp=viewport(x=x, y=y,
                            width=width, height=height,
                            name="face"),
        children=gList(rectGrob(vp="face"),
                       circleGrob(x=c(0.25, 0.75), 
                                  y=0.75, r=0.1, vp="face"),
                       linesGrob(x=c(0.33, 0.67), y=0.25,
                                 vp="face")))
}

faceD <- function(x, y, width, height) {
  grid.grabExpr({
                  pushViewport(viewport(x=x, y=y,
                                        width=width, 
                                        height=height))
                  grid.rect()
                  grid.circle(x=c(0.25, 0.75), 
                              y=0.75, r=0.1)
                  grid.lines(x=c(0.33, 0.67), y=0.25)
                  popViewport()
                })
}

drawDetails.face <- function(x, recording) {
  pushViewport(viewport(x=x$x, y=x$y,
                        width=x$width, height=x$height))
  grid.rect()
  grid.circle(x=c(0.25, 0.75), y=0.75, r=0.1)
  grid.lines(x=c(0.33, 0.67), y=0.25)
  popViewport()  
}

faceE <- function(x, y, width, height) {
  grob(x=x, y=y, width=width, height=height, cl="face")
}


figure8.25 <- function() {
faceA <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.rect()
  grid.circle(x=c(0.25, 0.75), y=0.75, r=0.1)
  grid.lines(x=c(0.33, 0.67), y=0.25)
  popViewport()
}

faceB <- function(x, y, width, height) {
  pushViewport(viewport(x=x, y=y, 
                        width=width, height=height))
  grid.draw(rectGrob())
  grid.draw(circleGrob(x=c(0.25, 0.75), y=0.75, r=0.1))
  grid.draw(linesGrob(x=c(0.33, 0.67), y=0.25))
  popViewport()
}


grid.newpage()
grid.draw(faceC(.5, .5, .5, .5))

grid.rect(gp=gpar(col="gray", fill=NA))



}
figure8.26 <- function() {
pdf("Figures/interactgrid-latticevps-%d.pdf", onefile=FALSE)
print(
xyplot(pressure ~ temperature, pressure)

)
showViewport(newpage=TRUE, 
             col="black",
             fill=rgb(.5, .5, .5, .2))
dev.off()
png("Web/interactgrid-latticevps%d.png")
print(
xyplot(pressure ~ temperature, pressure)

)
showViewport(newpage=TRUE, 
             col="black",
             fill=rgb(.5, .5, .5, .2))
dev.off()
system("rm Web/interactgrid-latticevps1.png")



}
figure8.27 <- function() {
pdf("Figures/interactgrid-latticeleaves-%d.pdf", onefile=FALSE)
print(
xyplot(pressure ~ temperature, pressure)

)
showViewport(newpage=TRUE, leaves=TRUE,
             col="black",
             fill=rgb(.5, .5, .5, .2))
dev.off()
png("Web/interactgrid-latticeleaves%d.png")
print(
xyplot(pressure ~ temperature, pressure)

)
showViewport(newpage=TRUE, leaves=TRUE,
             col="black",
             fill=rgb(.5, .5, .5, .2))
dev.off()
system("rm Web/interactgrid-latticeleaves1.png")




}
figure17.1 <- function() {

n <- 40
t <- seq(0, 2*pi, length=n)
x <- cos(t)
y <- sin(t)



for (i in 1:n) {
    plot.new()
    plot.window(c(-1, 1), c(-1, 1))
    lines(x, y)
    points(x[i], y[i], pch=16, cex=2)
    Sys.sleep(.05)
}



}
figure17.2 <- function() {





n <- 40
t <- seq(0, 2*pi, length=n)
x <- cos(t)
y <- sin(t)



orbit <- function() {
    par(pty="s", mar=rep(1, 4))
    for (i in 1:n) {
        plot.new()
        plot.window(c(-1, 1), c(-1, 1))
        lines(x, y)
        points(x[i], y[i], pch=16, cex=2)
    }
}



animation::ani.options(interval=0.05, outdir="orbitImages",
            filename="orbit.html")
animation::ani.start()
orbit()
animation::ani.stop()


}
figure17.3 <- function() {
plot(mpg ~ disp, mtcars)
points(mpg ~ disp, mtcars, subset=gear == 3, pch=16)



plot(mpg ~ disp, mtcars)
points(mpg ~ disp, mtcars, subset=gear == 4, pch=16)



}
figure17.4 <- function() {
trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
print(
xyplot(mpg ~ disp | factor(gear), mtcars, subset=gear != 5, pch=16)

)



}
figure17.5 <- function() {

mtcars$gear <- factor(mtcars$gear)
mtcars$cyl <- factor(mtcars$cyl)
gg <- rggobi::ggobi(mtcars)


}
figure17.6 <- function() {



}
figure17.7 <- function() {




gg <- rggobi::ggobi(mtcars)


}
figure17.8 <- function() {




iplots::iplot(mtcars$disp, mtcars$mpg)
iplots::ibar(mtcars$gear)


}
figure17.9 <- function() {




iplots::iplot(mtcars$disp, mtcars$mpg)
iplots::ibar(mtcars$gear)


iplots::iplot.set(1)



labels <- mapply("itext", 
                 mtcars$disp, mtcars$mpg, rownames(mtcars), 
                 MoreArgs=list(visible=FALSE), SIMPLIFY=FALSE)
olds <- NULL
while (!is.null(iplots::ievent.wait())) {
    if (iplots::iset.sel.changed()) {
        s <- iplots::iset.selected()
        if (length(s) > 1) {
            lapply(labels[s], iplots::iobj.opt, visible = TRUE)
        }
        if (length(olds) > 1) {
            lapply(labels[olds], iplots::iobj.opt, visible = FALSE)
        }
        olds <- s
    }
}


}
figure17.10 <- function() {



}
figure17.11 <- function() {



}
figure17.12 <- function() {
# library(latticist)
# latticist(mtcars, list(xvar="disp", yvar="mpg"), 
#           use.playwith=FALSE)
# 
# 
}
figure17.13 <- function() {

playwith::playwith(xyplot(mpg ~ disp, mtcars))
playwith::playwith(xyplot(qsec ~ wt, mtcars), 
         new=TRUE, link.to=playwith::playDevCur())


}
figure17.15 <- function() {
drawClock <- function(hour, minute) {
    t <- seq(0, 2*pi, length=13)[-13]
    x <- cos(t)
    y <- sin(t)

    grid.newpage()
    pushViewport(dataViewport(x, y, gp=gpar(lwd=4)))
    # Circle with ticks
    grid.circle(x=0, y=0, default.units="native", 
                r=unit(1, "native"))
    grid.segments(x, y, x*.9, y*.9, default.units="native")
    # Hour hand
    hourAngle <- pi/2 - (hour + minute/60)/12*2*pi
    grid.segments(0, 0, 
                  .6*cos(hourAngle), .6*sin(hourAngle), 
                  default.units="native", gp=gpar(lex=4))
    # Minute hand
    minuteAngle <- pi/2 - (minute)/60*2*pi
    grid.segments(0, 0, 
                  .8*cos(minuteAngle), .8*sin(minuteAngle), 
                  default.units="native", gp=gpar(lex=2))    
    grid.circle(0, 0, default.units="native", r=unit(1, "mm"),
                gp=gpar(fill="white"))
}






window <- gWidgets::gwindow("Clock")



allContent <- gWidgets::ggroup(container=window, horizontal=FALSE)



graphicTime <- gWidgets::ggraphics(container=allContent)



timeContent <- gWidgets::ggroup(container=allContent)



textLabel <- gWidgets::glabel("")



randomizeTime <- function(h, ...) {
    hour <- sample(1:12, 1)
    minute <- sample(seq(0, 55, 5), 1)
    drawClock(hour, minute)
    gWidgetsRGtk2::visible(textLabel) <- FALSE
    gWidgetsRGtk2::svalue(textLabel) <- paste(hour, 
                               sprintf("%02d", minute), 
                               sep=":")
}



reset <- gWidgets::gbutton("Randomize Time", 
                 handler=randomizeTime)



textButton <- gWidgets::gbutton("Show Time", 
                      handler=function(h, ...) {
                          gWidgetsRGtk2::visible(textLabel) <- TRUE
                      })


gWidgetsRGtk2::add(timeContent, reset)
gWidgetsRGtk2::add(timeContent, textButton)
gWidgetsRGtk2::add(timeContent, textLabel)



}
figure17.16 <- function() {




doc <- SVGAnnotation::svgPlot({ par(mfrow=c(2, 1), cex=.7,     
                     mar=c(5.1, 4.1, 1, 1))
                 plot(mpg ~ disp, mtcars, cex=2)
                 plot(qsec ~ wt, mtcars, cex=2) },
               width=4, height=8)
SVGAnnotation::linkPlots(doc)
XML::saveXML(doc, "linkedplots.svg")


}
figure1.1 <- function() {
plot(pressure)
text(150, 600, 
     "Pressure (mm Hg)\nversus\nTemperature (Celsius)")



}
figure14.1 <- function() {

colorado <- maptools::readShapeSpatial(system.file("extra", "10m-colorado.shp",
                                         package="RGraphics"))
par(mar=rep(0, 4))
sp::plot(colorado, col="gray")



}
figure14.2 <- function() {





par(mar=rep(0, 4))
maps::map(regions="Brazil", fill=TRUE, col="gray")

box("figure", col="gray", lwd=2)



}
figure14.3 <- function() {




brazil <- 
    maptools::readShapeSpatial(system.file("extra", "10m-brazil.shp",
                                 package="RGraphics"))



par(mar=rep(0, 4))
# Need to spec. sp:: here so that code works in 'RGraphics' package
sp::plot(brazil, col="gray")
box("figure", col="gray", lwd=2)



}
figure14.4 <- function() {




brazil <- 
    maptools::readShapeSpatial(system.file("extra", "10m-brazil.shp",
                                 package="RGraphics"))



print(
sp::spplot(brazil, "Regions", col.regions=gray(5:1/6))

)



}
figure14.5 <- function() {




brazil <- 
    maptools::readShapeSpatial(system.file("extra", "10m-brazil.shp",
                                 package="RGraphics"))



brazilRegions <- 
    maptools::readShapeSpatial(system.file("extra", 
                                 "10m_brazil_regions.shp",
                                 package="RGraphics"))


brazilCapitals <- 
    maptools::readShapeSpatial(system.file("extra",
                                 "10m_brazil_capitals.shp",
                                 package="RGraphics"))




print(
sp::spplot(brazil, "Regions", 
       col.regions=gray.colors(5, 0.8, 0.3),
       col="white", 
       panel=function(...) {
           sp::panel.polygonsplot(...)
           sp::sp.lines(brazilRegions, col="gray40")
           labels <- brazilCapitals$Name
           w <- stringWidth(labels)
           h <- stringHeight(labels)
           locs <- sp::coordinates(brazilCapitals)
           grid.rect(unit(locs[, 1], "native"),
                     unit(locs[, 2], "native"),
                     w, h, just=c("right", "top"),
                     gp=gpar(col=NA, fill=rgb(1, 1, 1, .5)))
           sp::sp.text(locs, labels, adj=c(1, 1))
           sp::sp.points(brazilCapitals, pch=21,
                     col="black", fill="white")
       })


)



}
figure14.6 <- function() {




marajo <- 
    maptools::readShapeSpatial(system.file("extra", "marajo.shp",
                                 package="RGraphics"))



par(mar=rep(0, 4))
sp::plot(marajo, col="gray", pbg="white")



}
figure14.7 <- function() {




iceland <- 
    maptools::readShapeSpatial(system.file("extra", "10m-iceland.shp",
                                 package="RGraphics"))



par(mar=rep(0, 4))
sp::plot(iceland, col="gray")



}
figure14.8 <- function() {




iceland <- 
    maptools::readShapeSpatial(system.file("extra", "10m-iceland.shp",
                                 package="RGraphics"))



sp::proj4string(iceland) <- sp::CRS("+proj=longlat +ellps=WGS84")



par(mar=rep(0, 4))
sp::plot(iceland, col="gray")



}
figure14.9 <- function() {




iceland <- 
    maptools::readShapeSpatial(system.file("extra", "10m-iceland.shp",
                                 package="RGraphics"))



sp::proj4string(iceland) <- sp::CRS("+proj=longlat +ellps=WGS84")



icelandMercator <- sp::spTransform(iceland, 
                               sp::CRS("+proj=merc +ellps=GRS80"))



par(mar=rep(0, 4))
sp::plot(iceland, col="gray80", border="white", lwd=3)
par(new=TRUE)
sp::plot(icelandMercator)



}
figure14.10 <- function() {








brazil <- 
    maptools::readShapeSpatial(system.file("extra", "10m-brazil.shp",
                                 package="RGraphics"))



sp::proj4string(brazil) <- sp::CRS("+proj=longlat +ellps=WGS84")



glines <- sp::gridlines(brazil)
glinesOrtho <- sp::spTransform(glines, sp::CRS("+proj=ortho"))
par(mar=rep(0, 4))
brazilOrtho <- sp::spTransform(brazil, sp::CRS("+proj=ortho"))
sp::plot(brazilOrtho, col="gray")
sp::plot(glinesOrtho, lty="dashed", add=TRUE)



}
figure14.11 <- function() {








brazil <- 
    maptools::readShapeSpatial(system.file("extra", "10m-brazil.shp",
                                 package="RGraphics"))



# Read in prepared raster
brazilRelief <- raster::raster(system.file("extra", "brazilRelief.tif",
                                   package="RGraphics"))



# Make PNG version for this one because otherwise it's TOO big
png("Figures/maps-brazilraster.png",
    width=900, height=900)
par(mar=rep(0, 4))
raster::image(brazilRelief, col=gray(0:255/255), maxpixels=1e6)
sp::plot(brazil, add=TRUE)
box(lwd=4)
dev.off()
system("cp Figures/maps-brazilraster.png Web/")


}
figure7.1 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.circle(name="circles", x=seq(0.1, 0.9, length=40), 
            y=0.5 + 0.4*sin(seq(0, 2*pi, length=40)),
            r=abs(0.1*cos(seq(0, 2*pi, length=40))))

grid.edit("circles", 
          gp=gpar(col=gray(c(1:20*0.04, 20:1*0.04))))

grid.remove("circles")




}
figure7.2 <- function() {
grid.rect(gp=gpar(col="gray"))
suffix <- c("even", "odd")
for (i in 1:8)
  grid.circle(name=paste("circle.", suffix[i %% 2 + 1], 
                         sep=""),
              r=(9 - i)/20, 
              gp=gpar(col=NA, fill=gray(i/10)))

grid.edit("circle.odd", gp=gpar(fill="gray10"), 
          global=TRUE)

grid.edit("circle", gp=gpar(col="gray", fill="gray90"), 
          grep=TRUE, global=TRUE) 




}
figure7.3 <- function() {
labels <- c("\"xaxis1\"\nxaxis gTree", "\"major\"\nlines grob", 
            "\"ticks\"\nlines grob", "\"labels\"\ntext grob")
names <- c("", "major", "ticks", "labels")
boxheight <- unit(2.5, "line")
boxwidth <- unit(1.2, "in")
pushViewport(viewport(layout=grid.layout(2, 3)))
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
grid.text(labels[1])
grid.lines(unit(0.5, "npc") + unit.c(-0.5*boxwidth, 0.5*boxwidth),
           0.5, gp=gpar(col="gray"))
grid.roundrect(height=boxheight, 
               width=boxwidth, # 1.2*stringWidth(labels[1]),
               r=unit(2, "mm"),
               gp=gpar(fill=NA))
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
grid.text(labels[2])
grid.lines(unit(0.5, "npc") + unit.c(-0.5*boxwidth, 0.5*boxwidth),
           0.5, gp=gpar(col="gray"))
grid.roundrect(height=boxheight, 
               width=boxwidth, # 1.2*stringWidth(labels[2]),
               r=unit(2, "mm"),
               gp=gpar(fill=NA))
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
grid.text(labels[3])
grid.lines(unit(0.5, "npc") + unit.c(-0.5*boxwidth, 0.5*boxwidth),
           0.5, gp=gpar(col="gray"))
grid.roundrect(height=boxheight, 
               width=boxwidth, # 1.2*stringWidth(labels[3]),
               r=unit(2, "mm"),
               gp=gpar(fill=NA))
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=3))
grid.text(labels[4])
grid.lines(unit(0.5, "npc") + unit.c(-0.5*boxwidth, 0.5*boxwidth),
           0.5, gp=gpar(col="gray"))
grid.roundrect(height=boxheight, 
               width=boxwidth, # 1.2*stringWidth(labels[4]),
               r=unit(2, "mm"),
               gp=gpar(fill=NA))
popViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
grid.move.to(x=0.5, y=unit(0.5, "npc") - 0.5*boxheight)
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
grid.line.to(x=0.5, y=unit(0.5, "npc") + 0.5*boxheight,
             arrow=arrow(angle=10, length=unit(3, "mm")))
popViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
grid.move.to(x=0.5, y=unit(0.5, "npc") - 0.5*boxheight)
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
grid.line.to(x=0.5, y=unit(0.5, "npc") + 0.5*boxheight,
             arrow=arrow(angle=10, length=unit(3, "mm")))
popViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
grid.move.to(x=0.5, y=unit(0.5, "npc") - 0.5*boxheight)
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=3))
grid.line.to(x=0.5, y=unit(0.5, "npc") + 0.5*boxheight,
             arrow=arrow(angle=10, length=unit(3, "mm")))
popViewport()



}
figure7.4 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(just="bottom", gp=gpar(cex=0.7)))
grid.xaxis(name="axis1", at=1:4/5)
grid.ls()

grid.edit("axis1", at=1:3/4)

grid.edit(gPath("axis1", "labels"), rot=45)

popViewport()



}
figure7.5 <- function() {
tg <- textGrob("sample text")
rg <- rectGrob(width=1.1*grobWidth(tg), 
               height=1.3*grobHeight(tg))
boxedText <- gTree(children=gList(tg, rg))



pushViewport(viewport(layout=grid.layout(1, 7,
                                         heights=unit(1.25, "in"),
                                         widths=unit(rep(c(1, 1.25), length=7),
                                                     rep(c("null", "in"),
                                                           length=7)))))
pushViewport(viewport(layout.pos.col=2, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
grid.draw(boxedText)

popViewport()
pushViewport(viewport(layout.pos.col=4, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
grid.draw(editGrob(boxedText, gp=gpar(col="gray")))

popViewport()
pushViewport(viewport(layout.pos.col=6, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
grid.draw(editGrob(boxedText, vp=viewport(angle=45),
                   gp=gpar(fontsize=18)))
              
popViewport()
popViewport()
                             


}
figure7.6 <- function() {
label <- textGrob("A\nPlot\nLabel ",
                  x=0, just="left")
x <- seq(0.1, 0.9, length=50)
y <- runif(50, 0.1, 0.9)
gplot <- 
  gTree(
    children=gList(rectGrob(gp=gpar(col="gray60",
                                    fill="white")),
                   linesGrob(x, y), 
                   pointsGrob(x, y, pch=16, 
                              size=unit(1.5, "mm"))),
    vp=viewport(width=unit(1, "npc") - unit(5, "mm"), 
                height=unit(1, "npc") - unit(5, "mm")))



layout <- grid.layout(1, 2,
                      widths=unit(c(1, 1), 
                                  c("null", "grobwidth"),
                                  list(NULL, label)))



grid.rect(gp=gpar(col="gray60", fill="gray90"))
pushViewport(viewport(layout=layout))
pushViewport(viewport(layout.pos.col=2))
grid.draw(label)
popViewport()
pushViewport(viewport(layout.pos.col=1))
grid.draw(gplot)
popViewport(2)




}
figure7.7 <- function() {
tg1 <- textGrob("Sample")
rg1 <- rectGrob(x=rep(0.5, 2),
                width=1.1*grobWidth(tg1), 
                height=1.3*grobHeight(tg1),
                gp=gpar(col=c("gray60", "white"), 
                        lwd=c(3, 1)))



pushViewport(viewport(layout=grid.layout(1, 7,
                                         heights=unit(1.25, "in"),
                                         widths=unit(rep(c(1, 1.25), length=7),
                                                     rep(c("null", "in"),
                                                           length=7)))))
pushViewport(viewport(layout.pos.col=2, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
grid.draw(tg1)
grid.draw(rg1)

popViewport()
pushViewport(viewport(layout.pos.col=4, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
pushViewport(viewport(gp=gpar(cex=2)))
grid.draw(tg1)
grid.draw(rg1)
popViewport()

popViewport()
pushViewport(viewport(layout.pos.col=6, gp=gpar(fill=NA)))
grid.rect(gp=gpar(col="gray", fill=NA))
pushViewport(viewport(gp=gpar(cex=2)))
grid.draw(tg1)
popViewport()
grid.draw(rg1)

popViewport()
popViewport()



}
figure7.8 <- function() {
tg1 <- textGrob("Sample", name="tg1")
rg1 <- rectGrob(width=1.1*grobWidth("tg1"), 
                height=1.3*grobHeight("tg1"),
                gp=gpar(col="gray60", lwd=3))
rg2 <- rectGrob(width=1.1*grobWidth(tg1), 
                height=1.3*grobHeight(tg1),
                gp=gpar(col="white"))



grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(cex=1.5, fill=NA)))
grid.draw(tg1)
grid.draw(rg1)
grid.draw(rg2)

grid.edit("tg1", grep=TRUE, global=TRUE, 
          label="Different text")

popViewport()



}
figure7.9 <- function() {
grid.rect(gp=gpar(col="gray"))
pushViewport(viewport(gp=gpar(fill=NA)))
grid.circle(.25, .5, r=unit(1, "mm"), 
            gp=gpar(fill="black"))
grid.text("A label", .75, .5)
grid.rect(.75, .5, 
          width=stringWidth("A label") + unit(2, "mm"),
          height=unit(1, "line"),
          name="labelbox")

grid.segments(.25, .5, 
              grobX("labelbox", 180), .5,
              arrow=arrow(angle=15, type="closed"),
              gp=gpar(fill="black"))





}
figure7.10 <- function() {
pushViewport(viewport(gp=gpar(fill=NA)))
vptop <- viewport(width=.9, height=.4, y=.75,
                  name="vptop")
vpbot <- viewport(width=.9, height=.4, y=.25,
                  name="vpbot")
pushViewport(vptop)
upViewport()
pushViewport(vpbot)
upViewport()

grid.rect(vp="vptop")
grid.lines(1:50/51, runif(50), vp="vptop")
grid.rect(vp="vpbot")
grid.lines(1:50/51, runif(50), vp="vpbot")

grid.null(x=.2, y=.95, vp="vptop", name="tl")
grid.null(x=.4, y=.95, vp="vptop", name="tr")
grid.null(x=.2, y=.05, vp="vpbot", name="bl")
grid.null(x=.4, y=.05, vp="vpbot", name="br")

grid.polygon(unit.c(grobX("tl", 0),
                    grobX("tr", 0),
                    grobX("br", 0),
                    grobX("bl", 0)),
             unit.c(grobY("tl", 0),
                    grobY("tr", 0),
                    grobY("br", 0),
                    grobY("bl", 0)),
             gp=gpar(col="gray", lwd=3))




}
figure7.11 <- function() {
grid.rect(gp=gpar(col="gray"))
grid.circle(r=0.3, gp=gpar(fill="gray80"), 
            name="mycircle")
grid.edit("mycircle", gp=gpar(lwd=5))
grid.edit("mycircle", gp=gpar(lty="dashed"))




}
figure7.12 <- function() {
angle <- seq(0, 2*pi, length=21)[-21]
x <- cos(angle)
y <- sin(angle)



trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
print(
xyplot(y ~ x, aspect=1, 
       xlab="displacement", 
       ylab="velocity")

)
grid.edit("[.]xlab$", grep=TRUE, 
          x=unit(1, "npc"), just="right",
          gp=gpar(fontfamily="mono"))
grid.edit("[.]ylab$", grep=TRUE, 
          y=unit(1, "npc"), just="right",
          gp=gpar(fontfamily="mono"))




}
figure7.13 <- function() {
mtcars2 <- mtcars
mtcars2$trans <- factor(mtcars$am, 
                        levels=0:1, 
                        labels=c("automatic", "manual"))
mtcars2$am <- NULL
mtcars2$vs <- NULL
mtcars2$drat <- NULL
mtcars2$carb <- NULL

# To keep R CMD check happy
mpg <- mtcars2$mpg



update_geom_defaults("smooth", aes(color="black"))
print(
ggplot(mtcars2, aes(x=disp, y=mpg)) +
    geom_point() +
    geom_smooth(method=lm)

)
downViewport("panel.3-4-3-4")
sline <- grid.get(gPath("smooth", "polyline"),
                  grep=TRUE)
grid.segments(.7, .8, 
              grobX(sline, 45), grobY(sline, 45),
              arrow=arrow(angle=10, type="closed"),
              gp=gpar(fill="black"))
grid.text("line of best fit", .71, .81,
          just=c("left", "bottom"))




}
figure1.14 <- function() {



grpkgs <- new("graphNEL", 
              nodes=c(
                # engine
                "grDevices", 
                # systems
                "graphics", "grid",
                # graphics-based packages
                "maps", "diagram", "plotrix", "gplots", "pixmap", 
                # grid-based packages
                "lattice", "ggplot2", "grImport", "gridBase", "vcd",
                # interface packages
                "rgl", "rggobi", "iplots",
                # devices
                "JavaGD", "Cairo", "tikzDevice"),
                # ggplot-based packages
                # "DescribeDisplay",
                # lattice-based packages
                # "latticist", "latticeExtra"),
              edgeL=list(
                grDevices=list(edges=c("graphics", "grid",
                                 "JavaGD", "Cairo", "tikzDevice")),
                graphics=list(edges=c("maps", "diagram", "plotrix",
                                "gplots", "pixmap", 
                                "gridBase")),
                grid=list(edges=c("lattice", "ggplot2", "grImport", "vcd",
                            "gridBase")),
                maps=list(),
                diagram=list(),
                plotrix=list(),
                gplots=list(),
                pixmap=list(),
                lattice=list(), # edges=c("latticist", "latticeExtra")), 
#                latticist=list(),
#                latticeExtra=list(),
                ggplot2=list(), # edges="DescribeDisplay"),
                # DescribeDisplay=list(),
                grImport=list(),
                gridBase=list(),
                vcd=list(),
                # Invisible links to tie interface packages together
                rgl=list(edges="rggobi"),
                rggobi=list(edges="iplots"),
                iplots=list(edges="rgl"),
#                GDD=list(edges=c("grDevices")),
#                JavaGD=list(edges=c("grDevices")),
#                Cairo=list(edges=c("grDevices")),
#                cairoDevice=list(edges=c("grDevices")),
#                tikzDevice=list(edges="grDevices")),
                JavaGD=list(),
                Cairo=list(),
                tikzDevice=list()),
              edgemode="directed")

# systemPkgs <- subGraph(c("graphics", "grid"), grpkgs)
# graphicsPkgs <- subGraph(c("maps", "diagram", "plotrix",
#                            "gplots", "pixmap"), grpkgs)
# gridPkgs <- subGraph(c("grid", "lattice", "ggplot2", "grImport"), grpkgs)
# 
# devicePkgs <- subGraph(c("GDD", "JavaGD", "Cairo", 
#                          "cairoDevice", "tikzDevice"),
#                        grpkgs)
interfacePkgs <- graph::subGraph(c("iplots", "rggobi", "rgl"), grpkgs)

ragraph <- Rgraphviz::agopen(grpkgs, name="whatever",
                  # layoutType="dot", 
                  layoutType="dot", 
                  # layoutType="twopi", 
                  attrs=list(
                    node=list(fontname="Helvetica", fontsize=10),
                    edge=list(arrowhead="none"),
                    # NOTE: size and margins controlled below in call to 'dot'
                    graph=list(
                      root="grDevices",
                      # ratio=3/4,
                      rankdir="LR")),
#                      compound=TRUE)),
                  subGList=list(
#                    list(graph=devicePkgs), 
#                    list(graph=systemPkgs),
#                    list(graph=graphicsPkgs),
#                    list(graph=gridPkgs),
                    list(graph=interfacePkgs)))

Rgraphviz::nodeDataDefaults(ragraph, "style") <- "filled"
Rgraphviz::nodeDataDefaults(ragraph, "fillcolor") <- "gray90"
Rgraphviz::nodeData(ragraph, c("grDevices", "graphics", "grid", "lattice", "ggplot2"), 
         "style") <- "filled"
Rgraphviz::nodeData(ragraph, c("grDevices", "graphics", "grid", "lattice", "ggplot2"), 
         "fillcolor") <- "gray70"

# clusterData(ragraph, 0, "label") <- "Devices"
# clusterData(ragraph, 1, "label") <- "Systems"
# clusterData(ragraph, 2, "label") <- "Graphics-based Packages"
# clusterData(ragraph, 3, "label") <- "Grid-based Packages"
Rgraphviz::clusterData(ragraph, 0, "style") <- "dashed"

# Edge from "grDevices" to "Interface Packages" cluster
# (Needs existing link from grDevices to rgl)
# edgeDataDefaults(ragraph, "lhead") <- NA
# edgeData(ragraph, "grDevices", "rgl", "lhead") <- "cluster_1" 

# Edges within "Interface Packages" cluster
Rgraphviz::edgeData(ragraph, "rgl", "rggobi", "style") <- "invis"
Rgraphviz::edgeData(ragraph, "rggobi", "iplots", "style") <- "invis"
Rgraphviz::edgeData(ragraph, "iplots", "rgl", "style") <- "invis"
# edgeData(ragraph, "pixmap", "rgl", "style") <- "invis"
# edgeData(ragraph, "gridBase", "rggobi", "style") <- "invis"
# edgeData(ragraph, "lattice", "iplots", "style") <- "invis"

Rgraphviz::toFile(ragraph, filename="grpkgs.dot", fileType="dot")
system("dot -Kneato grpkgs.dot -Tps -Gsize='8,8' -Gmargin=0 > organisation-graphicslevels.ps") 





}
figure10.1 <- function() {
names <- c("blank", "solid", "dashed", "dotted", "dotdash",
  "longdash", "twodash", "", "", "13", "F8", "431313", "22848222")
hgap <- unit(5, "mm")
pushViewport(viewport(layout=grid.layout(17, 5,
  widths=unit.c(unit(1, "strwidth", "Integer") + unit(3, "mm"), 
    hgap, unit(1.5, "inches"),
    hgap, max(unit(rep(1, length(names)), "strwidth", data=as.list(names))) + unit(3, "mm")),
  heights=unit(c(1.2, 1, 1, rep(1.2, 7), 1, 1, rep(1.2, 4), 1), "lines"))))
pushViewport(viewport(layout.pos.col=1:5, layout.pos.row=1:17))
grid.rect(width=1.3, height=1.3, gp=gpar(col="gray"))
popViewport()
for (i in 1:17) {
	if (i == 1) {
	  pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))  
	  grid.text("Integer", gp=gpar(fontface="bold"), just="right", x=1)
	  popViewport()
	  pushViewport(viewport(layout.pos.col=3, layout.pos.row=i))  
	  grid.text("Sample line", gp=gpar(fontface="bold"))
	  popViewport()
	  pushViewport(viewport(layout.pos.col=5, layout.pos.row=i))  
	  grid.text("String", gp=gpar(fontface="bold"), just="left", x=0)
	  popViewport()
	  pushViewport(viewport(layout.pos.row=i))  
	  grid.lines(c(0, 1), 0)
	  popViewport()
        } else if (i == 3) {
	    pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))  
	    grid.text("Predefined", just="left", x=0, 
              gp=gpar(fontface="italic"))
	    popViewport()          
        } else if (i == 12) {
	    pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))  
	    grid.text("Custom", just="left", x=0, 
              gp=gpar(fontface="italic"))
	    popViewport()          
	} else if ((i > 3 && i < 11) || 
                   (i > 12 && i < 17)) {
	  if (i < 11) {
	    pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))  
	    grid.text(i-4, just="right", x=1)
	    popViewport()
	  }
	  if (nchar(names[i-3])) {
		  pushViewport(viewport(layout.pos.col=5, layout.pos.row=i))  
		  grid.text(paste("\"", names[i-3], "\"", sep=""), x=0, just="left")
		  popViewport()
		  pushViewport(viewport(layout.pos.col=3, layout.pos.row=i))  
		  grid.lines(c(0, 1), 0.5, gp=gpar(lty=names[i-3], lwd=2))
		  popViewport()
	  }
        }
	if (i == 17) {
	    	  pushViewport(viewport(layout.pos.row=i))  
		  grid.lines(c(0, 1), 0)
		  popViewport()
	}
}
popViewport()



}
figure10.2 <- function() {
x <- c(.3, .7, .3)
y <- c(.2, .5, .8)
grid.rect(gp=gpar(col="gray"))
grid.lines(x, y, gp=gpar(lwd=40, lineend="square",
                   linejoin="mitre", col="black"))
grid.lines(x, y, gp=gpar(lwd=40, col="gray50"))
                   # lineend="round", linejoin="round"
grid.lines(x, y, gp=gpar(lwd=40, lineend="butt",
                   linejoin="bevel", col="gray80"))
grid.points(x, y, default.units="npc", pch=16, gp=gpar(cex=0.5))




}
figure10.3 <- function() {
ncol <- 6
nrow <- 5
grid.rect(gp=gpar(col="gray"))
for (i in 1:nrow) {
  for (j in 1:ncol) {
    x <- unit(j/(ncol+1), "npc")
    y <- unit(i/(nrow + 1), "npc")
    pch <- (i - 1)*ncol + j - 1
    if (pch > 25) 
      pch <- c("A", "b", ".", "#")[pch - 25]
    grid.points(x + unit(3, "mm"), y, 
      pch=pch, gp=gpar(fill="gray"))
    grid.text(pch, x - unit(3, "mm"), y, gp=gpar(col="gray"))
  }
}



}
figure10.4 <- function() {
h <- 0.01
drawexpr <- function(expr, y, exprFamily="CM") {
  grid.text(paste("expression(", expr, ")", sep=""), .5, y-h, 
            just="top", gp=gpar(fontfamily="mono", cex=0.75))
  grid.text(parse(text=expr), .5, y+h, 
            just="bottom", gp=gpar(fontfamily=exprFamily))
}
drawexpr("z[i] == sqrt(x[i]^2 + y[i]^2)", 1/5)
drawexpr("hat(beta) == (X^t * X)^{-1} * X^t * y", 2/5)
drawexpr("bar(x) == sum(frac(x[i], n), i==1, n)", 3/5)
drawexpr("paste(\"Temperature (\", degree, \"C) in 2003\")", 4/5,
         exprFamily="CM2")
grid.rect(gp=gpar(col="gray", fill=NA))



}
figure2.1 <- function() {
par(mfrow=c(2, 2), cex=0.6, mar=c(4, 4, 1, 1))
y <- rnorm(20)
plot(y, type="p")
plot(y, type="l")
plot(y, type="b")
plot(y, type="h")




}
figure2.2 <- function() {
  par(mfrow=c(2, 2), cex=0.6, mar=c(4, 4, 4, 2), mex=0.8)
  plot(lm.SR <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings),
       id.n=1, cex.caption=0.8, which=1:4,
       panel=function(...) { panel.smooth(..., col.smooth="gray") })



}
figure2.3 <- function() {




subset <- sample(1:150, 20)
cS <- as.character(Sp <- iris$Species[subset])
cS[Sp == "setosa"] <- "S"
cS[Sp == "versicolor"] <- "V"
cS[Sp == "virginica"] <- "g"
ai <- cluster::agnes(iris[subset, 1:4])

par(mfrow=c(2, 1), cex=0.5, pty="s")
plot(ai, which=1, col=c("gray90", "gray"), labels = cS)
plot(ai, which=2, labels = cS)



}
figure2.4 <- function() {
plottitle <- function(plotfun, funarg, outer=FALSE, cex=.7, line=1) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    mtext(paste(plotfun, "(", 
                paste(rep(" ", nca), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
plot2title <- function(plotfun, funarg, 
                       extrafn, extraarg, 
                       outer=FALSE, cex=.7, line=.5) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    ncep <- nchar(extrafn)
    ncea <- nchar(extraarg)
    mtext(paste(plotfun, 
                "(",  paste(rep(" ", nca), collapse=""),
                ")\n", 
                extrafn, "(",
                paste(rep(" ", ncea), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " \n", 
                paste(rep(" ", ncep + 1), collapse=""),
                extraarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
dohplot <- function(plotfn, ..., funarg, 
                    extrafn=NULL, extraarg=NULL, namefudge=FALSE,
                    main="", xlab="", ylab="", axes=FALSE, box=TRUE) {
    if (is.null(xlab) || is.null(ylab)) {
        do.call(plotfn,
                list(..., main=""))
    } else if (is.null(axes)) {
        do.call(plotfn,
                list(..., main="", xlab="", ylab=""))
    } else {
        do.call(plotfn,
                list(..., main="", axes=FALSE, xlab="", ylab=""))
    }
    if (is.null(extrafn)) {
        plottitle(plotfn, funarg)
    } else {
        plot2title(if (namefudge) paste(" ", plotfn, sep="") else plotfn, 
                   funarg, extrafn, extraarg)
    }
    if (box)
        box() # col="gray")
}



par(mfrow=c(3, 4), mar=c(1, 1, 3, 1), mex=.7, mgp=c(3, 100, 100))
dohplot("plot", (1:10)^2, funarg="numeric",
        xlim=c(0, 11), ylim=c(-10, 110))
dohplot("plot", table(rep(1:3, 1:3)), funarg="table", 
        lwd=2, xlim=c(0, 4), ylim=c(0, 4))
# Empty
plot.new()
plot.new()
dohplot("barplot", table(rep(1:3, 1:3)), funarg="", 
        extrafn="plot", extraarg="factor",
        xlim=c(-1, 5), ylim=c(0, 4), names.arg="")
dohplot("pie", c(1, 2, 4), funarg="", col=gray(1:3/4), cex=.7,
        labels="", axes=NULL)
dohplot("dotchart", 1:3, 
        funarg="numeric", pch=21, bg="gray",  
        lcolor="black", xlim=c(0, 4))
# Empty
plot.new()
dohplot("boxplot", (1:10)^2, funarg="numeric", 
        col="gray", ylim=c(-10, 110))
dohplot("hist", (1:100)^2, funarg="", col="gray", 
        breaks=6,
        xlim=c(-1000, 11000), ylim=c(0, 50))
dohplot("stripchart", (1:10)^2, funarg="numeric", 
        method="stack",
        cex=1, xlim=c(-10, 110), ylim=c(-1, 3), pch=21, bg="gray")
# stem()
plot.new()
txt <- capture.output(stem((1:10)^2))[-2]
text(.05, (1:length(txt))/(length(txt) + 1), txt, adj=0, family="mono", cex=.7)
box() # col="gray")
plottitle("stem", "")



}
figure2.5 <- function() {
plottitle <- function(plotfun, funarg, outer=FALSE, cex=.7, line=1) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    mtext(paste(plotfun, "(", 
                paste(rep(" ", nca), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
plot2title <- function(plotfun, funarg, 
                       extrafn, extraarg, 
                       outer=FALSE, cex=.7, line=.5) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    ncep <- nchar(extrafn)
    ncea <- nchar(extraarg)
    mtext(paste(plotfun, 
                "(",  paste(rep(" ", nca), collapse=""),
                ")\n", 
                extrafn, "(",
                paste(rep(" ", ncea), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " \n", 
                paste(rep(" ", ncep + 1), collapse=""),
                extraarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
dohplot <- function(plotfn, ..., funarg, 
                    extrafn=NULL, extraarg=NULL, namefudge=FALSE,
                    main="", xlab="", ylab="", axes=FALSE, box=TRUE) {
    if (is.null(xlab) || is.null(ylab)) {
        do.call(plotfn,
                list(..., main=""))
    } else if (is.null(axes)) {
        do.call(plotfn,
                list(..., main="", xlab="", ylab=""))
    } else {
        do.call(plotfn,
                list(..., main="", axes=FALSE, xlab="", ylab=""))
    }
    if (is.null(extrafn)) {
        plottitle(plotfn, funarg)
    } else {
        plot2title(if (namefudge) paste(" ", plotfn, sep="") else plotfn, 
                   funarg, extrafn, extraarg)
    }
    if (box)
        box() # col="gray")
}



set.seed(1500)
# mgp draws the axes miles off the page
par(mfrow=c(4, 4), mar=c(1, 1, 3, 1), mex=.7, mgp=c(3, 100, 100))
dohplot("plot", 1:10, (1:10)^2, funarg="num,num", 
        pch=21, bg="gray", 
        xlim=c(0, 11), ylim=c(-10, 110))
x <- rnorm(10000)
dohplot("smoothScatter", x, x + rnorm(10000)/3,  funarg="",
        nbin=64, colramp=function(n) { gray(n:1/(n + 1)) },
        xlim=c(-5, 5), ylim=c(-5, 5))
x <- sample(1:4, 20, replace=TRUE)
y <- x + sample(0:1, 20, replace=TRUE)
dohplot("sunflowerplot", x, y,
        funarg="", seg.col="black", size=.07,
        xlim=c(0, 5), ylim=c(0, 6))
# Empty gap 
plot.new()
dohplot("boxplot", list((1:10)^2, 120 - (1:10)^2), funarg="list", 
        extrafn="plot", extraarg="fac,num", col="gray", boxwex=0.5,
        ylim=c(-10, 130))
dohplot("barplot", rbind(1:3, (1:3)^2), funarg="matrix",
        xlim=c(0, 4), ylim=c(0, 13))
dohplot("barplot", rbind(1:3, (1:3)^2), funarg="matrix",
        beside=TRUE,
        xlim=c(0, 10), ylim=c(0, 11))
# Empty gap for dotchart
plot.new()
fig <- par("fig")
dohplot("stripchart", list((1:10)^2, 140 - (1:10)^2), funarg="list",
        extrafn="plot", extraarg="num,fac",
        xlim=c(-10, 150), ylim=c(0, 3), pch=21, bg="gray", cex=1)
dohplot("spineplot", 
        rep(1:3, each=6), 
        factor(c(rep(1:3, 3:1), rep(1:3, 2), rep(1:3, 1:3))),
        funarg="num,fac", box=FALSE)
dohplot("cdplot", 
        rep(1:3, each=6), 
        factor(c(rep(1:3, 3:1), rep(1:3, 2), rep(1:3, 1:3))),
        funarg="", box=FALSE)
# Empty gap 
plot.new()
dohplot("spineplot", 
        factor(rep(1:3, each=6)), 
        factor(c(rep(1:3, 3:1), rep(1:3, 2), rep(1:3, 1:3))),
        funarg="fac,fac", off=5,
        extrafn="plot", extraarg="fac,fac",
        namefudge=TRUE,
        box=FALSE)
dohplot("assocplot", 
        table(rep(1:2, each=3),
              c(rep(1:2, 1:2), rep(1:2, 2:1))),
        funarg="",
        col=c("black", "gray"), axes=NULL)
dohplot("fourfoldplot", 
        table(rep(1:2, each=3),
              c(rep(1:2, 1:2), rep(1:2, 2:1))),
        color=gray(1:2/3),
        # NOTE: can't make 'space' too small or font size of labels
        # goes to zero and get error from ghostscript
        funarg="", xlab=NULL, box=FALSE, space=0.03)
dohplot("mosaicplot", 
        table(factor(rep(1:3, each=6)), 
              factor(c(rep(1:3, 3:1), rep(1:3, 2), rep(1:3, 1:3)))),
        funarg="", off=10,
        extrafn="plot", extraarg="table",
        cex.axis=.1, box=FALSE)
# Put dotchart in gap
par(fig=c(fig[1] - .1, fig[2:4]), new=TRUE)
dotdata <- rbind(1:3, (1:3)^2) # rbind(table(gy), table(gx))
colnames(dotdata) <- rep("", 3)
dohplot("dotchart", dotdata, funarg="matrix",
        labels="", pch=c(16, 21), bg="gray",
        lcolor="black", xlim=c(0, 13), box=FALSE)



}
figure2.6 <- function() {
plottitle <- function(plotfun, funarg, outer=FALSE, cex=.7, line=1) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    mtext(paste(plotfun, "(", 
                paste(rep(" ", nca), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
plot2title <- function(plotfun, funarg, 
                       extrafn, extraarg, 
                       outer=FALSE, cex=.7, line=.5) {
    ncp <- nchar(plotfun)
    nca <- nchar(funarg)
    ncep <- nchar(extrafn)
    ncea <- nchar(extraarg)
    mtext(paste(plotfun, 
                "(",  paste(rep(" ", nca), collapse=""),
                ")\n", 
                extrafn, "(",
                paste(rep(" ", ncea), collapse=""),
                ")", sep=""),
          family="mono", cex=cex, line=line, font=2, outer=outer)
    mtext(paste(paste(rep(" ", ncp + 1), collapse=""),
                funarg, " \n", 
                paste(rep(" ", ncep + 1), collapse=""),
                extraarg, " ", sep=""),
          family="mono", col="gray60", cex=cex, line=line, font=2, outer=outer)
}
dohplot <- function(plotfn, ..., funarg, 
                    extrafn=NULL, extraarg=NULL, namefudge=FALSE,
                    main="", xlab="", ylab="", axes=FALSE, box=TRUE) {
    if (is.null(xlab) || is.null(ylab)) {
        do.call(plotfn,
                list(..., main=""))
    } else if (is.null(axes)) {
        do.call(plotfn,
                list(..., main="", xlab="", ylab=""))
    } else {
        do.call(plotfn,
                list(..., main="", axes=FALSE, xlab="", ylab=""))
    }
    if (is.null(extrafn)) {
        plottitle(plotfn, funarg)
    } else {
        plot2title(if (namefudge) paste(" ", plotfn, sep="") else plotfn, 
                   funarg, extrafn, extraarg)
    }
    if (box)
        box() # col="gray")
}



# mgp draws the axes miles off the page
par(mfrow=c(3, 4), mar=c(1, 1, 3, 1), mex=.7, mgp=c(3, 100, 100))
mdf <- cbind(3:6, (3:6)^2, (3:6)^3)
names(mdf) <- c("Y1", "Y2", "Y3")
aaa <- seq(0, pi, length=10)
xxx <- rep(aaa, 10)
yyy <- rep(aaa, each=10)
zzz <- sin(xxx) + sin(yyy)
# gap for pairs(matrix)
plot.new()
dohplot("matplot", mdf[order(mdf[, 1]), ], funarg="", 
        pch=21:23, bg=c("white", "black", "gray"), type="o",
        col="black", xlim=c(0, 6), ylim=c(-10, 230))
df <- USJudgeRatings[1:6, ]
rownames(df) <- NULL
dohplot("stars", df, funarg="", ncol=2, lwd=1,
        len=.8, col.stars=rep("gray", 13), mar=c(1, 1, 3, 1))
# gap
plot.new()
dohplot("image", matrix(zzz, ncol=10), funarg="", col=gray(1:12/13))
dohplot("contour", matrix(zzz, ncol=10), funarg="", 
        levels=seq(0, 2, .25), labcex=.4)
# gap for filled.contour(matrix)
plot.new()
dohplot("persp", matrix(zzz, ncol=10), funarg="",
        theta=30, phi=45, col="gray")
dohplot("symbols", xxx, yyy, funarg="",
        circles=zzz, inches=.03, axes=NULL)
# gap for coplot(y ~ x | g)
plot.new()
dohplot("mosaicplot", 
        table(factor(rep(1:3, each=6)), 
              factor(c(rep(1:3, 3:1), rep(1:3, 2), rep(1:3, 1:3)))),
        funarg="", cex.axis=.1, off=10,
        box=FALSE)



}
figure2.7 <- function() {
par(mfrow=c(2, 2), mar=c(2.5, 2, 1, 1), cex=0.6)
boxplot(decrease ~ treatment, data = OrchardSprays,
        log = "y", col="light gray")
boxplot(decrease ~ treatment, data = OrchardSprays,
        log = "y", col="light gray", 
        boxwex=0.5)

par(las=2, xpd=NA)
barplot(VADeaths[1:2,], angle = c(45, 135), 
        density = 20, col = "gray",
        names=c("RM", "RF", "UM", "UF"))
barplot(VADeaths[1:2,], angle = c(45, 135), 
        density = 20, col = "gray",
        names=c("RM", "RF", "UM", "UF"),
        horiz=TRUE)




}
figure2.8 <- function() {
par(mfrow=c(2, 2), mar=c(2, 2, 1, 1), cex=0.6)
y <- rnorm(20)
plot(y, type="l", lwd=3)
plot(y, type="l", col="gray")
plot(y, type="l", lty="dashed")
plot(y, type="l", ylim=c(-4, 4))




}
figure2.9 <- function() {
par(cex=.5)
plot(function(x) { 
         sin(x)/x 
     }, 
     from=-10*pi, to=10*pi, 
     xlab="", ylab="", n=500)




par(mfrow=c(1, 2))
par(mar=c(7, 0, 3, 1))
par(mex=0.7)

hc <- hclust(dist(USArrests), "ave")
dend1 <- as.dendrogram(hc)
dend2 <- cut(dend1, h=70)
par(cex=0.7)
par(mar=c(1, 0, 2, 8))
#  dend2$lower is *NOT* a dendrogram, but a list of .. :
plot(dend2$lower[[3]], 
  horiz = TRUE, type = "tr", axes=FALSE, cex=0.8)
par(mar=c(8, 0, 2, 0))
# "inner" and "leaf" edges in different type & color :
plot(dend2$lower[[2]], 
     edgePar = list(col=c("black", "gray")), edge.root=TRUE, 
     axes=FALSE, cex=0.8)



}
figure4.1 <- function() {

trellis.par.set(list(dot.symbol=list(pch=1)))
print(
xyplot(pressure ~ temperature, pressure)

)



}
figure4.2 <- function() {
tplot <- xyplot(pressure ~ temperature, pressure)




trellis.par.set(list(dot.symbol=list(pch=1)))
print(
xyplot(pressure ~ temperature, pressure,
       type="o", pch=16, lty="dashed", 
       main="Vapor Pressure of Mercury")


)



}
figure4.3 <- function() {
x <- 1:5
y <- 1:5
g <- factor(1:5)
types <- c("barchart", "bwplot", "densityplot", "dotplot",
           "histogram", "qqmath", "stripplot", "qq",
           "xyplot", "levelplot", "contourplot",
           "cloud", "wireframe", "splom", "parallel")
angle <- seq(0, 2*pi, length=19)[-19]
xx <- cos(angle)
yy <- sin(angle)
gg <- factor(rep(1:3, each=6))

aaa <- seq(0, pi, length=10)
xxx <- rep(aaa, 10)
yyy <- rep(aaa, each=10)
zzz <- sin(xxx) + sin(yyy)


doplot <- function(name, ...) {
  do.call(name, 
          list(..., scales=list(draw=FALSE), xlab=NULL, ylab=NULL,
               strip=function(which.panel, ...) { 
                       grid.rect(gp=gpar(fill="gray90")); grid.text(name) 
                     }))
}
plot <- vector("list", 15)
plot[[1]] <- doplot("barchart", y ~ g | 1)
plot[[2]] <- doplot("bwplot", yy ~ gg | 1, 
                    par.settings=list(box.umbrella=list(lwd=0.5)))
plot[[3]] <- doplot("densityplot", ~ yy | 1)
plot[[4]] <- doplot("dotplot", g ~ y | 1)
plot[[5]] <- doplot("histogram", ~ xx | 1)
plot[[6]] <- doplot("qqmath", ~ yy | 1)
plot[[7]] <- doplot("stripplot", yy ~ gg | 1)
plot[[8]] <- doplot("qq", gg ~ yy | rep(1, 18), subset=gg != 3)
plot[[9]] <- doplot("xyplot", xx ~ yy | 1)
plot[[10]] <- doplot("levelplot", zzz ~ xxx + yyy | 1, colorkey=FALSE)
plot[[11]] <- doplot("contourplot", zzz ~ xxx + yyy | 1, labels=FALSE, cuts=8)
plot[[12]] <- doplot("cloud", zzz ~ xxx + yyy | 1, zlab=NULL, zoom=0.9, 
                     par.settings=list(box.3d=list(lwd=0.1)))
plot[[13]] <- doplot("wireframe", zzz ~ xxx + yyy | 1, zlab=NULL, zoom=0.9,
                     drape=TRUE, par.settings=list(box.3d=list(lwd=0.1)),
                     colorkey=FALSE)
plot[[14]] <- doplot("splom", ~ data.frame(x=xx[1:10], y=yy[1:10]) | 1, 
                     pscales=0)
plot[[15]] <- doplot("parallel", ~ as.data.frame(split(yy, gg)) | 1)

grid.newpage()
pushViewport(viewport(layout=grid.layout(4, 4)))
for (i in 1:15) {
  pushViewport(viewport(layout.pos.col=((i - 1) %% 4) + 1,
                        layout.pos.row=((i - 1) %/% 4) + 1))
  print(plot[[i]], newpage=FALSE, 
        panel.width=list(1.025, "inches"),
        panel.height=list(1.025, "inches"))
  popViewport()
}
popViewport()
 


}
figure4.4 <- function() {

print(
xyplot(mpg ~ disp, data=mtcars)

)



}
figure4.5 <- function() {

trellis.par.set(list(dot.symbol=list(pch=1)))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0, ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars)

)



}
figure4.6 <- function() {

trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0, ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
xyplot(mpg ~ disp, data=mtcars,
       group=gear, 
       auto.key=list(space="right"),
       par.settings=list(superpose.symbol=list(pch=c(1, 3, 16),
                           fill="white")))

)



}
figure4.7 <- function() {

trellis.par.set(list(dot.symbol=list(pch=1)))
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars,
       layout=c(1, 3), aspect=1)

)



}
figure4.8 <- function() {

trellis.par.set(list(fontsize=list(text=10)))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0, ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
plot1 <- xyplot(mpg ~ disp, data=mtcars, 
                aspect=1, xlim=c(65, 480), ylim=c(9, 35),
                subset=gear == 5)
plot2 <- xyplot(mpg ~ disp, data=mtcars, 
                aspect=1, xlim=c(65, 480), ylim=c(9, 35),
                subset=gear == 4)
plot3 <- xyplot(mpg ~ disp, data=mtcars, 
                aspect=1, xlim=c(65, 480), ylim=c(9, 35),
                subset=gear == 3)
print(plot1, position=c(0, 2/3, 1, 1), more=TRUE)
print(plot2, position=c(0, 1/3, 1, 2/3), more=TRUE)
print(plot3, position=c(0, 0, 1, 1/3))




}
figure4.9 <- function() {

trellis.par.set(list(fontsize=list(text=10)))
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars,
       layout=c(3, 1), aspect=1,
       scales=list(y=list(at=seq(10, 30, 10))),
       ylab="miles per gallon",
       xlab=expression(paste("displacement (", inch^3, ")")))

)



}
figure4.10 <- function() {

trellis.par.set(list(fontsize=list(text=10)))
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars,
       layout=c(3, 1), aspect=1,
       panel=function(...) {
           panel.xyplot(...)
           panel.abline(h=29, lty="dashed")
           panel.text(470, 29.5, "efficiency criterion",
                      adj=c(1, 0), cex=.7)
       })

)



}
figure4.11 <- function() {

trellis.par.set(list(fontsize=list(text=10)))
gray.colors <- function(n) { 
    adjustcolor(gray(n:1/n), alpha.f=.7) 
}
print(
xyplot(mpg ~ disp | factor(gear), data=mtcars,
       layout=c(3, 1), aspect=1,
       panel=function(x, y, ...) {
           panel.lmline(x, y)
           panel.xyplot(x, y, ...)
       })

)



}
figure4.12 <- function() {

trellis.par.set(list(fontsize=list(text=9, points=8)))
show.settings()



}
figure4.13 <- function() {
doplot <- function(name, ...) {
  do.call(name, 
          list(..., scales=list(draw=FALSE), xlab=NULL, ylab="",
               strip=function(which.panel, ...) { 
                       grid.rect(gp=gpar(fill="gray90")); grid.text(name) 
                     }))
}

plot <- vector("list", 4)
plot[[1]] <- doplot("ecdfplot", ~ rnorm(10) | 1)
plot[[2]] <- doplot("rootogram", ~ rpois(50, 10) | 1, 
                    dfun=function(x) dpois(x, 10))
plot[[3]] <- doplot("segplot", factor(1:5) ~ rnorm(5) + rnorm(5) | 1)
plot[[4]] <- doplot("tileplot", 
                    1:10 ~ rnorm(10) + rnorm(10) | 1,
                    colorkey=FALSE, aspect="fill")


grid.newpage()
pushViewport(viewport(layout=grid.layout(1, 4)))
for (i in 1:4) {
  pushViewport(viewport(layout.pos.col=i,
                        layout.pos.row=1))
  print(plot[[i]], newpage=FALSE, 
        panel.width=list(1.025, "inches"),
        panel.height=list(1.025, "inches"))
  popViewport()
}
popViewport()



}
figure12.1 <- function() {
TitanicDF <- as.data.frame(Titanic)
TitanicList <- lapply(TitanicDF[1:4], rep, TitanicDF$Freq)
TitanicSets <- 
    data.frame(passenger=TitanicList$Class != "Crew",
               adult=TitanicList$Age == "Adult",
               male=TitanicList$Sex == "Male",
               survivor=TitanicList$Survived == "Yes")
head(TitanicSets)











gplots::venn(TitanicSets[1:2])



par(mar=rep(2, 4))
plot(venneuler::venneuler(TitanicSets[1:2]), 
     col=hcl(0, 0, c(60, 80), .5),
     alpha=NA, border="black")




gplots::venn(TitanicSets[1:3])



par(mar=rep(2, 4))
plot(venneuler::venneuler(TitanicSets[1:3]), 
     col=hcl(0, 0, seq(40, 80, 20), .5),
     alpha=NA, border="black")




pdf("Figures/special-venn-3-%d.pdf", onefile=FALSE,
    width=6, height=6)
gplots::venn(TitanicSets)

dev.off()
png("Figures/special-venn-3-%d.png", width=240, height=240, pointsize=8)
gplots::venn(TitanicSets)

dev.off()


par(mar=rep(2, 4))
plot(venneuler::venneuler(TitanicSets[1:4]), 
     col=hcl(0, 0, seq(20, 80, 20), .5),
     alpha=NA, border="black")




}
figure12.2 <- function() {








TeachingDemos::faces(USJudgeRatings[1:5, ], nrow=1, ncol=5)



TeachingDemos::faces2(USJudgeRatings[1:5, ], nrows=1, ncols=5, scale="all")



symbols::symbol(USJudgeRatings[1:25, ], type="face")



}
figure12.3 <- function() {

data("soil", envir=environment())
place2region <- as.data.frame(rbind(c("Cnt1", "Lima"),
                                    c("Cnt2", "Lima"),
                                    c("Cnt3", "Lima"),
                                    c("Chz", "Ancash"),
                                    c("Chmar", "Huanuco"),
                                    c("Hco1", "Huanuco"),
                                    c("Hco2", "Huanuco"),
                                    c("Hco3", "Huanuco"),
                                    c("Hyo1", "Junin"),
                                    c("Hyo2", "Junin"),
                                    c("Namora", "Cajamarca"),
                                    c("SR1", "Junin"),
                                    c("SR2", "Junin")))
soils <- merge(soil, place2region,
               by.x="place", by.y="V1")[c("sand", "slime", "clay")]
names(soils) <- c("sand", "silt", "clay")















vcd::ternaryplot(soils, col="black", 
            grid_color="black", labels_color="black")
 



plotrix::triax.plot(soils,  cex.ticks=.5)



TTsoils <- soils
names(TTsoils) <- c("SAND", "SILT", "CLAY")
soiltexture::TT.plot(tri.data=TTsoils)



}
figure12.4 <- function() {
hourSpeed <- aggregate(RGraphics::hourlySpeed["Speed"], 
                       list(hour=RGraphics::hourlySpeed$hour),
                       mean)
head(hourSpeed)











trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
print(
with(RGraphics::wind9am,
     polarFreq(data.frame(ws=Speed, wd=Dir, date=Date),
               cols=gray(10:1/11), border.col="black"))

      )



plotrix::polar.plot(hourSpeed$Speed, hourSpeed$hour * 15,
           start=90, clockwise=TRUE, lwd=5,
           label.pos=seq(15, 360, 15), labels=1:24,
           radial.lim=c(0, 4.5))


}
figure12.5 <- function() {
 



station22254dir <- with(RGraphics::wind9am, Dir[Station == 22254])



station22254 <- circular::circular(station22254dir, 
                         units="degrees",
                         zero=pi/2, rotation="clock")



windHours <- circular::circular(RGraphics::hourlySpeed$hour,
                      units="hours", 
                      zero=pi/2, rotation="clock")



par(mar=rep(2, 4), xpd=NA)
plot(station22254, stack=TRUE, sep=.06)




par(mar=rep(2, 4), xpd=NA)
plot(density(station22254, bw=45), 
     main="", xlab="", ylab="")




par(mar=rep(1, 4), xpd=NA)
circular::rose.diag(station22254, bins=36, prop=3)




par(mar=rep(1, 4), xpd=NA)
plot(windHours, col=NA, shrink=1.2, axes=FALSE)
lines(windHours, 
      0.5*RGraphics::hourlySpeed$Speed/max(RGraphics::hourlySpeed$Speed),
      nosort=TRUE, lty="dotted", join=FALSE)
circular::axis.circular(template="clock24")




}
figure12.6 <- function() {





print(
with(RGraphics::wind9am,
     windRose(data.frame(ws=Speed, wd=Dir, 
                         date=Date, station=factor(Station)),

              paddle=FALSE, type="station", width=2))
)


}
figure12.7 <- function() {




data("NHANES", envir=environment())
plot(Serum.Iron ~ Transferin, NHANES)



trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
print(
hexbin::hexbinplot(Serum.Iron ~ Transferin, NHANES)

)



trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
print(
hexbin::hexbinplot(Serum.Iron ~ Transferin | Sex, NHANES)

)



}
figure16.1 <- function() {

tetra <- function() {
t1 <- rgl::tetrahedron3d()
t2vb <- t1$vb
t2vb[1, ] <- -3
t2 <- rgl::tmesh3d(t2vb, t1$it)
plane <- rgl::qmesh3d(rbind(rep(-3.01, 4),
                       c(-2, -2, 2, 2),
                       c(-3, 3, 3, -3),
                       rep(1, 4)),
                 matrix(1:4, ncol=1))
rgl::open3d(windowRect=c(0, 0, 600, 600))
# clear3d()
rgl::shade3d(plane, color="white", specular="black")
rgl::wire3d(plane)
rgl::wire3d(t1, lwd=3)
rgl::wire3d(t2, lwd=3)
rgl::segments3d(rbind(t2$vb[1, t2$it], 
                 t1$vb[1, t1$it]),
           rbind(t2$vb[2, t2$it], 
                 t1$vb[2, t1$it]),
           rbind(t2$vb[3, t2$it], 
                 t1$vb[3, t1$it]),
           col="gray", lwd=3)
rgl::view3d(40, -30)
}


t1 <- rgl::cube3d()
t1tube <- t1
t1tube$ib <- t1tube$ib[, -(3:4)]
t2vb <- t1$vb
t2vb[1, ] <- -5
t2 <- rgl::qmesh3d(t2vb, t1$ib)
plane <- rgl::qmesh3d(rbind(rep(-5.01, 4),
                       c(-2, -2, 2, 2),
                       c(-3, 3, 3, -3),
                       rep(1, 4)),
                 matrix(1:4, ncol=1))
rgl::open3d(windowRect=c(0, 0, 600, 600))
# clear3d()
rgl::shade3d(plane, color="white", 
        ambient="white", specular="white", emission="white")
rgl::wire3d(plane)
rgl::shade3d(t1tube, color="white", specular="black")
rgl::wire3d(t1, lwd=3)
rgl::wire3d(t2, lwd=3)
rgl::segments3d(rbind(t2$vb[1, t2$ib[, 4]], 
                 t1$vb[1, t1$ib[, 4]]),
           rbind(t2$vb[2, t2$ib[, 4]], 
                 t1$vb[2, t1$ib[, 4]]),
           rbind(t2$vb[3, t2$ib[, 4]], 
                 t1$vb[3, t1$ib[, 4]]),
           col="gray", lwd=3)
rgl::view3d(55, -20, fov=0)
rgl::rgl.postscript("Figures/threed-3dproj.eps")
system("epstopdf Figures/threed-3dproj.eps")
system("convert Figures/threed-3dproj.pdf Web/threed-3dproj.png")



}
figure16.2 <- function() {

tetra <- function() {
t2vb <- t1$vb
t2vb[1, ] <- -3
t2vb[2, c(1, 4)] <- t2vb[2, c(1, 4)]*.8
t2vb[3, c(1, 4)] <- t2vb[3, c(1, 4)]*.8
t2vb[2, 2:3] <- t2vb[2, 2:3]*.6
t2vb[3, 2:3] <- t2vb[3, 2:3]*.6
t2 <- rgl::tmesh3d(t2vb, t1$it)
t3vb <- t1$vb
t3vb[1, ] <- -10
t3vb[2, ] <- 0
t3vb[3, ] <- 0
t3 <- rgl::tmesh3d(t3vb, t1$it)
rgl::open3d(windowRect=c(0, 0, 600, 600))
# clear3d()
rgl::shade3d(plane, color="white", specular="black")
rgl::wire3d(plane)
rgl::wire3d(t1, lwd=3)
rgl::wire3d(t2, lwd=3)
rgl::segments3d(rbind(t3$vb[1, t3$it], 
                 t1$vb[1, t1$it]),
           rbind(t3$vb[2, t3$it], 
                 t1$vb[2, t1$it]),
           rbind(t3$vb[3, t3$it], 
                 t1$vb[3, t1$it]),
           col="gray", lwd=3)
rgl::shade3d(rgl::translate3d(rgl::scale3d(rgl::cube3d(), .1, .1, .1), -10, 0, 0))
rgl::view3d(50, -20, zoom=.8)
}

t1 <- rgl::cube3d()
t1tube <- t1
t1tube$ib <- t1tube$ib[, -(3:4)]
t2vb <- t1$vb
t2vb[2, t1$ib[, 4]] <- t2vb[2, t1$ib[, 4]]*.4
t2vb[3, t1$ib[, 4]] <- t2vb[3, t1$ib[, 4]]*.4
t2vb[2, t1$ib[, 3]] <- t2vb[2, t1$ib[, 3]]*.6
t2vb[3, t1$ib[, 3]] <- t2vb[3, t1$ib[, 3]]*.6
t2vb[1, ] <- -4.5
t2 <- rgl::qmesh3d(t2vb, t1$ib)
t2tube <- t2
t2tube$ib <- t2tube$ib[, -(3:4)]
t3 <- rgl::translate3d(rgl::scale3d(rgl::cube3d(), .1, .1, .1), -10, 0, 0)
plane <- rgl::qmesh3d(rbind(rep(-4.51, 4),
                       c(-2, -2, 2, 2),
                       c(-3, 3, 3, -3),
                       rep(1, 4)),
                 matrix(1:4, ncol=1))
rgl::open3d(windowRect=c(0, 0, 600, 300))
# clear3d()
rgl::shade3d(plane, color="white", 
        ambient="white", specular="white", emission="white")
rgl::wire3d(plane)
rgl::shade3d(t1tube, color="white", specular="black")
rgl::wire3d(t1, lwd=3)
rgl::shade3d(t2tube, color="white", specular="black")
rgl::wire3d(t2, lwd=3)
rgl::segments3d(rbind(-10, 
                 t1$vb[1, t1$ib[, 4]]),
           rbind(0, 
                 t1$vb[2, t1$ib[, 4]]),
           rbind(0,
                 t1$vb[3, t1$ib[, 4]]),
           col="gray", lwd=3)
rgl::shade3d(t3)
rgl::view3d(50, -15, fov=0, zoom=.7)
rgl::rgl.postscript("Figures/threed-3dvp.eps")
system("epstopdf Figures/threed-3dvp.eps")
system("convert Figures/threed-3dvp.pdf Web/threed-3dvp.png")



}
figure16.3 <- function() {
quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



par(mar=rep(0, 4))
persp(quakeDens)




}
figure16.4 <- function() {
quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



par(mar=rep(0, 4))
persp(quakeDens, scale=FALSE, expand=0.02,
      theta=60, d=.1, r=.1,
      xlab="longitude", ylab="latitude", zlab="")




}
figure16.5 <- function() {
quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



par(mar=rep(0, 4))
zinterp <- with(quakeDens,
                z[-1, -1] + z[-1, -ncol(z)] + 
                z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
persp(quakeDens, scale=FALSE, expand=0.02,
      theta=60, d=.1, r=.1, axes=FALSE, box=FALSE,
      col=gray(.4 + 1:6/10)[cut(zinterp, 6)])




}
figure16.6 <- function() {
quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



shallowCantyQuakes <- subset(cantyQuakes, DEPTH < 20)



trellis.device("pdf", color=FALSE,
               file="Figures/threed-cloud%d.pdf", onefile=FALSE)
for (i in seq(40, 80, 20)) {
    print(cloud(-DEPTH ~ LONG + LAT, shallowCantyQuakes,
                xlim=c(171, 173), ylim=c(-44.5, -42.5),
                pch=16, col=rgb(0, 0, 0, .5),
                screen=list(z=i, x=-70)))
}

dev.off()



}
figure16.7 <- function() {
quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



shallowCantyQuakes <- subset(cantyQuakes, DEPTH < 20)







par(lab=c(3, 3, 0))
s3d <- with(shallowCantyQuakes,
            scatterplot3d(-DEPTH ~ LONG + LAT,
                          angle=30, scale.y=0.45, type="n",
                          pch=16, color=rgb(0, 0, 0, .5),
                          x.ticklabs=pretty(LONG, 3),
                          grid=FALSE, zlim=c(-20, 0)))

quakeDensXY <- MASS::kde2d(shallowCantyQuakes$LONG, 
                     shallowCantyQuakes$LAT, n=30)
lapply(contourLines(quakeDensXY, nlevels=8),
       function(cl) {
           polygon(s3d$xyz.convert(cl$x, cl$y, 
                                   rep(-20, length(cl$x))),
                   lwd=.5, col=gray(.8 - cl$level/20),
                   border=NA)
       })

quakeDensXZ <- MASS::kde2d(shallowCantyQuakes$LONG, 
                     -shallowCantyQuakes$DEPTH, n=30)
lapply(contourLines(quakeDensXZ, nlevels=8),
       function(cl) {
           polygon(s3d$xyz.convert(cl$x, 
                                   rep(-43.2, length(cl$x)),
                                   cl$y),
                   lwd=.5, col=gray(.8 - 2*cl$level),
                   border=NA)
       })
quakeDensYZ <- MASS::kde2d(shallowCantyQuakes$LAT, 
                     -shallowCantyQuakes$DEPTH, n=30)
lapply(contourLines(quakeDensYZ, nlevels=8),
       function(cl) {
           polygon(s3d$xyz.convert(rep(171.5, length(cl$x)),
                                 cl$x, cl$y),
                   lwd=.5, col=gray(.8 - cl$level),
                   border=NA)
       })
lapply(contourLines(quakeDensYZ, nlevels=8),
       function(cl) {
           polygon(s3d$xyz.convert(rep(173.5, length(cl$x)),
                                   cl$x, cl$y),
                   lwd=.5, col=gray(.8 - cl$level),
                   border=NA)
       })

with(shallowCantyQuakes,
     s3d$points3d(-DEPTH ~ LONG + LAT, pch=16,
                  col=rgb(0, 0, 0, .3)))

s3d$box()



}
figure16.8 <- function() {




quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



rgl::open3d(windowRect=c(0, 0, 900, 450))
# clear3d("all")
rgl::persp3d(quakeDens$x, quakeDens$y, quakeDens$z, 
        aspect=c(1, 0.55, .2), col="white", box=FALSE,
        axes=FALSE, xlab="", ylab="", zlab="")

rgl::par3d(userMatrix=rgl::rotationMatrix(-80/180*pi, 1, 0, 0)%*%
                 rgl::rotationMatrix(-65/180*pi, 0, 0, 1),
      zoom=.5)
rgl::snapshot3d("Figures/threed-rglpersp.png")
system("cp Figures/threed-rglpersp.png Web/")



}
figure16.9 <- function() {








quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



shallowCantyQuakes <- subset(cantyQuakes, DEPTH < 20)



d <- with(shallowCantyQuakes, 
          {
              misc3d::kde3d(LONG, LAT, -DEPTH, 
                    h=c(.1, .1, 2), n = 30)
          })



rgl::open3d(windowRect=c(0, 0, 900, 900))
with(shallowCantyQuakes, 
     {
         rgl::plot3d(LONG, LAT, -DEPTH, 
                aspect=c(1, 0.55, 1), 
                axes=TRUE, box=FALSE,
                xlab="", ylab="", zlab="")
         misc3d::contour3d(d$d, c(.4, .1), d$x, d$y, d$z,
                   color=rep("gray80", 2), 
                   color2="gray", specular="black",
                   engine="rgl", add=TRUE, alpha=.5)
     })

rgl::par3d(userMatrix=rgl::rotationMatrix(-60/180*pi, 1, 0, 0)%*%
                 rgl::rotationMatrix(-40/180*pi, 0, 0, 1),
      zoom=.9)
rgl::snapshot3d("Figures/contour3d.png")
system("cp Figures/contour3d.png Web/threed-contour3d.png")


}
figure16.10 <- function() {




rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d("all")
rgl::light3d()
rgl::material3d(shininess=100, specular="black")
# Head
radius <- function(d) {
   pchisq(d^2, 3)
}
rgl::shade3d(rgl::ellipse3d(diag(3), level=radius(1),
                   centre=c(0, 0, 1)),
         color="yellow")
# Neck
# logo is 100x76
png("rlogoExtended.png",
    width=500, height=250)
grid.rect(gp=gpar(col=NA, fill="yellow"))

rlogo <- as.raster(png::readPNG(system.file("extra", "Rlogo.png", 
                                       package="RGraphics")))
rlogo[rlogo == "#FFFFFF"] <- "yellow"
grid.raster(rlogo, x=.6, y=.01, width=.08, just=c("bottom"))
dev.off()
rgl::shade3d(rgl::cylinder3d(cbind(0, 0, c(-1.4, 1)),
                   e1=cbind(0, 0, 1),
                   e2=cbind(1, 0, 0),
                   sides=100),
        color="yellow",
        texture="rlogoExtended.png",
        texcoords=list(x=rep(seq(1, 0, length.out=101), each=4)[-c(1:2, 403:404)],
                       y=rep(c(0, 1, 1, 0), 50)))
old <- function() {
rgl::shade3d(rgl::cylinder3d(cbind(0, 0, c(-1.3, 1)),
                   e1=cbind(0, 0, 1),
                   e2=cbind(1, 0, 0),
                   sides=100),
        color="yellow")
}
# Eyes
eyeball <- rgl::ellipse3d(diag(3), level=radius(.4))
rgl::shade3d(rgl::translate3d(eyeball, .8, .35, .7),
         color="white")
rgl::shade3d(rgl::translate3d(eyeball, .8, -.35, .7),
         color="white")
# Translate radius of eye, rotate, translate position of eye
pupil <- rgl::rotate3d(rgl::translate3d(rgl::ellipse3d(diag(3),
                                        level=radius(.05)),
                              .4, 0, 0),
                  30/180*pi, 0, 0, 1)
rgl::shade3d(rgl::translate3d(pupil, .8, .35, .7),
         color="black")
rgl::shade3d(rgl::translate3d(pupil, .8, -.35, .7),
         color="black")
# points3d(1.21, c(-.35, .35), .7, cex=3)
# Nose
rgl::shade3d(rgl::cylinder3d(cbind(c(1, 1.3), 0, .3),
                    radius=.2,
                    e1=cbind(1, 0, 0),
                    e2=cbind(0, 1, 0),
                    sides=100),
         color="yellow")
rgl::shade3d(rgl::ellipse3d(diag(3), level=radius(.2),
                   centre=c(1.3, 0, .3)),
         color="yellow")
# Mouth
rgl::shade3d(rgl::ellipse3d(diag(3), level=radius(.8),
                   centre=c(.6, 0, -.5)),
         color="tan")
angle <- seq(-65, 65, length=30)/180*pi
rgl::lines3d(.6 + .81*cos(angle), .81*sin(angle), -.5, lwd=3)
# Hair on top
angle <- seq(15, 165, length=30)/180*pi
rgl::lines3d(.2, .7*cos(angle), 1.5 + .7*sin(angle), lwd=3)
rgl::lines3d(-.2, .7*cos(angle), 1.5 + .7*sin(angle), lwd=3)
# Hair on sides
rgl::lines3d(seq(.5, -.5, length=5), -1, rep(c(.3, .8), length=5), lwd=3)
rgl::lines3d(seq(.5, -.5, length=5), 1, rep(c(.3, .8), length=5), lwd=3)

rgl::par3d(userMatrix=rgl::rotationMatrix(-pi/2, 1, 0, 0)%*%
      rgl::rotationMatrix(-50/180*pi, 0, 0, 1)%*%
      rgl::rotationMatrix(10/180*pi, 1, 0, 0))

rgl::snapshot3d("homer.png")



}
figure16.11 <- function() {




quakes <- read.csv("Quake/earthquakes.csv")
NZquakes <- quakes[c("LAT", "LONG", "MAG", "DEPTH")]



cantyQuakes <- quakes[quakes$LAT < -42.4 & quakes$LAT > -44 & 
               quakes$LONG > 171 & quakes$LONG < 173.5, ]

quakeDens <- MASS::kde2d(cantyQuakes$LONG, cantyQuakes$LAT, n=30)



shallowCantyQuakes <- subset(cantyQuakes, DEPTH < 20)



with(shallowCantyQuakes,
     vrmlgen::cloud3d(LONG, LAT, -DEPTH,
             filename="vrmlgen.wrl",
             cols="white"))




}
figure13.1 <- function() {
xmm <- read.table(file.path("XMM-Newton", "XMM-Newton.txt"),
                  header=TRUE)



counts <- sort(table(RGraphics::xmm$Category))

par(mfrow=c(1, 2), mar=c(3, 3, 2, 2))

barplot(counts)
dotchart(counts)




print(barchart(counts), pos=c(0, 0, .5, 1),
      more=TRUE)
print(dotplot(counts), pos=c(.5, 0, 1, 1))



grid.newpage()
catSort <- data.frame(Category=factor(RGraphics::xmm$Category, levels=names(counts)))
pushViewport(viewport(x=0, width=.5, just="left"))
print(ggplot(catSort) +
      geom_bar(aes(x=Category)), newpage=FALSE)
popViewport()
pushViewport(viewport(x=1, width=.5, just="right"))
catCounts <- data.frame(Category=factor(names(counts), levels=names(counts)),
                        Count=counts)
print(ggplot(catCounts) +
      geom_point(aes(y=Category, x=Count)),
      newpage=FALSE)
popViewport()




}
figure13.2 <- function() {




vcd::spine(Priority ~ Duration, RGraphics::xmm)



durn <- RGraphics::xmm$Duration/1000
vcd::cd_plot(Priority ~ durn, RGraphics::xmm, xlab="Duration (1000s)")



}
figure13.3 <- function() {
trellis.par.set(theme = canonical.theme("postscript", color=FALSE))
catTab <- table(RGraphics::xmm$Schedule, RGraphics::xmm$Priority)
print(barchart(prop.table(catTab, margin=1), col=gray(1:3/4)),
      pos=c(0, 0, .5, 1), more=TRUE)
print(barchart(prop.table(catTab, margin=1), col=gray(1:3/4), stack=FALSE),
      pos=c(.5, 0, 1, 1))



grid.newpage()
pushViewport(viewport(x=0, width=.5, just="left"))
print(
ggplot(as.data.frame(prop.table(catTab, margin=1))) +
    geom_bar(aes(x=Var1, fill=Var2, y=Freq), 
             stat="identity", col="black") +
    scale_fill_manual(values=gray(1:3/3)),
      newpage=FALSE)
popViewport()
pushViewport(viewport(x=1, width=.5, just="right"))
print(
ggplot(as.data.frame(prop.table(catTab, margin=1))) +
    geom_bar(aes(x=Var1, fill=Var2, y=Freq), 
             stat="identity", col="black", position="dodge") +
    scale_fill_manual(values=gray(1:3/3)),
      newpage=FALSE)
popViewport()




}
figure13.4 <- function() {




vcd::mosaic(Priority ~ Schedule, RGraphics::xmm)



vcd::mosaic(nObs ~ Schedule + Priority, RGraphics::xmm,
       labeling_args=list(rot_labels=c(right=0), 
         offset_labels=c(right=-.5),
         just_labels=c(right="left")),
       margin=c(right=4))




}
figure13.5 <- function() {




grid.rect(gp=gpar(col=NA, fill="gray"))
vcd::tile(nObs ~ Schedule + Priority, RGraphics::xmm,
     tile_type="area",
     shade=TRUE, 
     gp=gpar(lwd=2, fill="white"), 
     pos_labels=c(left="left", top="left", right="left"), 
     just_labels=c(left="left", top="left", right="left"),
     pop=FALSE, newpage=FALSE)

downViewport("cell:Schedule=fixed,Priority=C,nObs=multiple")
grid.circle(0, 0, r=unit(1, "mm"))
upViewport(0)

downViewport("cell:Schedule=fixed,Priority=C,nObs=single")
grid.circle(0, 0, r=unit(1, "mm"))
upViewport(0)
downViewport("cell:Schedule=free,Priority=A,nObs=multiple")
grid.circle(0, 0, r=unit(1, "mm"))
upViewport(0)
downViewport("cell:Schedule=free,Priority=B,nObs=multiple")
grid.circle(0, 0, r=unit(1, "mm"))
upViewport(0)




vcd::doubledecker(nObs ~ Schedule + Priority, RGraphics::xmm,
             dep_varname=FALSE,
             gp=gpar(fill=c("gray90", "gray")),
             offset_labels=c(right=-.5),
             margins=c(bottom=3, left=1, top=1, right=5))




}
figure13.6 <- function() {
xmm <- read.table(file.path("XMM-Newton", "XMM-Newton.txt"),
                  header=TRUE)







pairs(vcd::structable(nObs ~ Priority + Schedule, RGraphics::xmm),
      space=.15)



}
figure13.7 <- function() {
vcd::cotabplot(~ Schedule + Priority | nObs, RGraphics::xmm)



}
figure13.8 <- function() {




# Need the .1 cos the handling of zero cells seems off
# Also need the custom shading to produce gray-scale
vcdExtra::mosaic3d(vcd::structable(~ Priority + Schedule + nObs, RGraphics::xmm) + .1, 
         shading=function(x) { "gray" })



}
