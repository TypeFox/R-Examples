`plot3dstatic.anacor` <-
function(x, plot.type = "jointplot", plot.dim = c(1,2,3), col.r = "RED", col.c = "BLUE", arrows = TRUE, 
                                 main, xlab, ylab, zlab, xlim, ylim, zlim, ...)
{
#produces static 3D-scatterplot
#plot.type can be of "jointplot", "colplot", "rowplot".

options(locatorBell = FALSE)
if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !")
if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
pd1 <- plot.dim[1]
pd2 <- plot.dim[2]
pd3 <- plot.dim[3]
if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")

if (missing(xlab)) xlab <- paste("Dimension",pd1)
if (missing(ylab)) ylab <- paste("Dimension",pd2)
if (missing(zlab)) zlab <- paste("Dimension",pd3)

#--------------------------------- rowplot ----------------------------------

if (plot.type == "rowplot") {
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(main)) main1 <- paste("Row Plot") else main1 <- main

  #if (missing(col.r)) col.r <- 1

  xa <- x$row.scores[,plot.dim]
  if (missing(xlim)) xlim <- range(xa[,1])
  if (missing(ylim)) ylim <- range(xa[,2])
  if (missing(zlim)) zlim <- range(xa[,3])

  pr <- scatterplot3d(xa[,1], xa[,2], xa[,3], type = "n", main = main1, xlab = xlab, ylab = ylab,
                      zlab = zlab, xlim = xlim, ylim = ylim, zlim = zlim,...)
  text(pr$xyz.convert(xa[,1], xa[,2], xa[,3]), labels = rownames(x$row.scores), col = col.r)
  if (arrows) for (i in 1:nrow(xa)) pr$points3d(rbind(xa[i,],c(0,0,0)), col = col.r, type="l", lty=1)
}
#---------------------------------- end rowplot --------------------------------

#--------------------------------- colplot ----------------------------------

if (plot.type == "colplot") {
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(main)) main1 <- paste("Column Plot") else main1 <- main

  #if (missing(col.c)) col.c <- 2

  xa <- x$col.scores[,plot.dim]
  if (missing(xlim)) xlim <- range(xa[,1])
  if (missing(ylim)) ylim <- range(xa[,2])
  if (missing(zlim)) zlim <- range(xa[,3])

  pr <- scatterplot3d(xa[,1], xa[,2], xa[,3], type = "n", main = main1, xlab = xlab, ylab = ylab,
                      zlab = zlab, xlim = xlim, ylim = ylim, zlim = zlim,...)
  text(pr$xyz.convert(xa[,1], xa[,2], xa[,3]), labels = rownames(x$col.scores), col = col.c)
  if (arrows) for (i in 1:nrow(xa)) pr$points3d(rbind(xa[i,],c(0,0,0)), col = col.c, type="l", lty=1)
}
#---------------------------------- end colplot --------------------------------

#--------------------------------- jointplot ----------------------------------

if (plot.type == "jointplot") {
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(main)) main1 <- paste("Joint Plot") else main1 <- main

  #if (missing(col.c)) col.c <- 2
  #if (missing(col.r)) col.r <- 1
  

  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
  
  if (missing(xlim)) xlim <- range(c(xa[,1],ya[,1]))
  if (missing(ylim)) ylim <- range(c(xa[,2],ya[,2]))
  if (missing(zlim)) zlim <- range(c(xa[,3],ya[,3]))

  xy <- rbind(xa,ya)
  pr <- scatterplot3d(xy, type = "n", main = main1, xlab = xlab, ylab = ylab,
                      zlab = zlab, xlim = xlim, ylim = ylim, zlim = zlim,...)
  text(pr$xyz.convert(xa), labels = rownames(xa), col = col.r)
  text(pr$xyz.convert(ya), labels = rownames(ya), col = col.c)
  
  if (arrows) {
    for (i in 1:nrow(xa)) pr$points3d(rbind(xa[i,],c(0,0,0)), col = col.r, type="l", lty=1)
    for (i in 1:nrow(ya)) pr$points3d(rbind(ya[i,],c(0,0,0)), col = col.c, type="l", lty=1)
  }
}
#---------------------------------- end jointplot --------------------------------

}
