`plot3d.anacor` <-
function(x, plot.type = "jointplot", plot.dim = c(1,2,3), col.r = "RED", col.c = "BLUE", arrows = TRUE,
                          xlab, ylab, zlab, main, ...)
{


# x         ... object of class "anacor"
# plot.type ... string with either "rowplot", "colplot", "jointplot"
# plot3d.dim ... dimensions to be plotted (vector of length 3),
#               vector of length 2 for "rowplot", "colplot" and "jointplot", ignored for "benzplot".
# wline     ... line thickness for jointplot (for all others ignored)

bg.image <- "particle"
n <- dim(x$tab)[1]
m <- dim(x$tab)[2]
if (x$ndim < 3) stop("No 3D plot can be produced for ndim < 3!")

if (plot.type == "rowplot") {
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(main)) main <- paste("Row plot for",x$datname)
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(zlab)) zlab <- paste("Dimension",plot.dim[3])
  #if (missing(col.r)) col.r <- 2

  xa <- x$row.scores[,plot.dim]
  
  path <- paste("textures/",bg.image,".png",sep="")
  rgl.open()
  rgl.bg(sphere=TRUE, texture=system.file(path, package="rgl"),back = "filled", ...)
  text3d(xa, texts = rownames(xa), col = col.r, ...)
  axes3d(c('x','y','z'), ...)
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, ...)
  if (arrows) for (i in 1:nrow(x$row.scores)) rgl.lines(rbind(c(0,0,0),x$row.scores[i,]), col = col.r) 
}
#---------------------------------- end rowplot --------------------------------


#----------------------------------- colplot -----------------------------------
if (plot.type == "colplot") {
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(main)) main1 <- paste("Column plot for",x$datname)
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(zlab)) zlab <- paste("Dimension",plot.dim[3])
  #if (missing(col.c)) col.c <- 2

  ya <- x$col.scores[,plot.dim]
  path <- paste("textures/",bg.image,".png",sep="")
  rgl.open()
  rgl.bg(sphere=TRUE, texture=system.file(path, package="rgl"),back = "filled", ...)
  text3d(ya,texts = rownames(ya), col = col.c, ...)
  axes3d(c('x','y','z'), ...)
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, ...)
  if (arrows) for (i in 1:nrow(x$col.scores)) rgl.lines(rbind(c(0,0,0),x$col.scores[i,]), col = col.c)
}
#------------------------------- end colplot -----------------------------------


#--------------------------------- jointplot -----------------------------------
#jointplot makes a labelled biplot, using one of four different scalings.
#graph plot is optionally added

if (plot.type == "jointplot") {
  wlines <- 1
  if (length(plot.dim) != 3) stop("plot.dim must be of length 3")
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(zlab)) zlab <- paste("Dimension",plot.dim[3])
  if (missing(main)) main <- paste("Joint plot for",x$datname)
  #if (missing(col.r)) col.r <- 1
  #if (missing(col.c)) col.c <- 2

  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
  rg <- round(wlines*x$tab/max(x$tab))
  path <- paste("textures/",bg.image,".png",sep="")
  
  rgl.open()
  rgl.bg(sphere=TRUE, texture=system.file(path, package="rgl"),back = "filled")
  text3d(xa,texts = rownames(xa), col = col.r)
  text3d(ya,texts = rownames(ya), col = col.c)
  axes3d(c('x','y','z'))
  title3d(xlab = xlab, ylab = ylab, zlab = zlab)
  if (arrows) {
    for (i in 1:nrow(x$row.scores)) rgl.lines(rbind(c(0,0,0),x$row.scores[i,]), col = col.r)
    for (i in 1:nrow(x$col.scores)) rgl.lines(rbind(c(0,0,0),x$col.scores[i,]), col = col.c)
  }
}

#--------------------------------- end jointplot -------------------------------
}

