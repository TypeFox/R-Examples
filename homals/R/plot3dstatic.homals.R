plot3dstatic.homals <- function(x, plot.dim = c(1,2,3), plot.type = "jointplot", var.subset, 
                                 main, type, xlab, ylab, zlab, xlim, ylim, zlim, ...)
{
#produces static 3D-scatterplot
#plot.type can be of "objplot", "catplot", "labplot", "starplot", "jointplot".

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

nvar <- dim(x$dframe)[2]
if (missing(var.subset)) var.subset <- 1:nvar

x1 <- x$objscores[,pd1]
y1 <- x$objscores[,pd2]
z1 <- x$objscores[,pd3]

#----------------------------- object plot -------------------------------------
if (plot.type == "objplot") {

  if (missing(xlim)) xlim <- range(x1)
  if (missing(ylim)) ylim <- range(y1)
  if (missing(zlim)) zlim <- range(z1)

  if (missing(main)) main1 <- "Object score plot" else main1 <- main

  pr <- scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab, ylab = ylab,
                      zlab = zlab, xlim = xlim, ylim = ylim, zlim = zlim, ...)
  text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$objscores))
}
#----------------------------- end object plot ---------------------------------

#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- t(sapply(x$loadings, function(xy) xy[1,c(pd1,pd2,pd3)]))  #first solution only
  if (missing(main)) main1 <- "Loadings plot" else main1 <- main
 
  pr <- scatterplot3d(xycoor,  main = main1, xlab = xlab, ylab = ylab, zlab = zlab, type = "n",...)
  pr$points3d(xycoor, col = "RED")

  for (i in 1:nvar) 
    pr$points3d(rbind(xycoor[i,c(pd1,pd2,pd3)],c(0,0,0)), col = "BLUE", type="l", lty=1)
  identify(pr$xyz.convert(xycoor), labels = rownames(xycoor), col = "RED")
    
}
#-------------------------------- end loadplot ---------------------------------


#--------------------------------------labplot----------------------------------
#plot labeled object scores (for each variable separately)

if (plot.type == "labplot") {

  if (missing(xlim)) xlim <- range(x1)
  if (missing(ylim)) ylim <- range(y1)
  if (missing(zlim)) zlim <- range(z1)
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Labplot for",colnames(x$dframe[i]))  else main1 <- main
    par("ask" = TRUE)
    pr <- scatterplot3d(x1, y1, z1, type = "n", xlim = xlim, ylim = ylim, zlim = zlim, 
                  xlab = xlab, ylab = ylab, zlab = zlab, main = main1, ...)
    text(pr$xyz.convert(x1, y1, z1), labels = as.vector(x$dframe[,i]))
  }   
} 
#-----------------------------------end labplot---------------------------------

#----------------------------------catplot--------------------------------------
#plots the rank-restricted category quantifications for each variable

if (plot.type == "catplot") {

  if (missing(type)) type <- "b"
  if (missing(xlim)) xlim <- range(sapply(x$catscores, function(zz) range(zz[,pd1])))
  if (missing(ylim)) ylim <- range(sapply(x$catscores, function(zz) range(zz[,pd2])))
  if (missing(ylim)) zlim <- range(sapply(x$catscores, function(zz) range(zz[,pd3])))

  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Category plot for",colnames(x$dframe[i]))  else main1 <- main
    par("ask"=TRUE)
    pr <- scatterplot3d(x$catscores[[i]][,c(pd1,pd2,pd3)], type = type, xlim = xlim, ylim = ylim,
                        main = main1, xlab = xlab, ylab = ylab, zlab = zlab,  
                        col.grid = "white", ...)
    text(pr$xyz.convert(x$catscores[[i]][,c(pd1,pd2,pd3)]), labels = levels(x$dframe[,i]), 
         pos = 3)
    pr$plane3d(c(0,0,0),col = "gray", lty.box = "solid")
  }
}
#----------------------------------end catplot----------------------------------

if (plot.type == "starplot") {
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Star plot for", colnames(x$dframe[i])) else main1 <- main
    
    par("ask" = TRUE)
    pr <- scatterplot3d(x1, y1, z1, color = "RED", pch = 24, main = main1, 
                        xlab = xlab, ylab = ylab, zlab = zlab, type = "n", ...)    
    y <- computeY(x$dframe[,i], x$objscores[,c(pd1,pd2,pd3)])
    pr$points3d(y, col = "RED", pch = 24)
    text(pr$xyz.convert(y), labels = levels(x$dframe[,i]), pos = 3, col="RED")
    pr$points3d(x1, y1, z1, col = "BLUE")
    #text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$dframe), pos = 3, col="BLUE")

    for (j in 1:length(x$dframe[,i])) 
      pr$points3d(rbind(x$objscores[j,c(pd1,pd2,pd3)],y[x$dframe[,i][j],]), col = "BLUE", type="l", lty=1)
    identify(pr$xyz.convert(x1, y1, z1), labels = rownames(x$dframe), col = "BLUE") 
  }
}
#----------------------------------end starplot---------------------------------

if (plot.type == "jointplot") {
  xylist <- lapply(x$catscores, apply, 2, range)         
  xytab <- sapply(xylist, function(yy) yy[,c(pd1,pd2,pd3)])
  xmin <- min(xytab[1,], x$objscores[,pd1])
  xmax <- max(xytab[2,], x$objscores[,pd1])
  ymin <- min(xytab[3,], x$objscores[,pd2])
  ymax <- max(xytab[4,], x$objscores[,pd2])
  zmin <- min(xytab[5,], x$objscores[,pd3])
  zmax <- max(xytab[6,], x$objscores[,pd3])
    
  if (missing(xlim)) xlim <- c(xmin, xmax)
  if (missing(ylim)) ylim <- c(ymin, ymax)
  if (missing(zlim)) zlim <- c(zmin, zmax)
  
  if (missing(main)) main <- "Joint Plot"
  pr <- scatterplot3d(x1, y1, z1, type = "n", main = main, xlab = xlab, ylab = ylab, 
                      zlab = zlab, xlim = xlim, ylim = ylim, zlim = zlim, ...)           #draw scores
  text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$dframe), col = 1)
  catcol <- rainbow(ncol(x$dframe))
  
  catleg <- NULL
  for (j in var.subset)
  {
    text(pr$xyz.convert(x$catscores[[j]][,c(pd1,pd2,pd3)]), labels = rownames(x$catscores[[j]]), col = catcol[j])
    catleg <- c(catleg, catcol[j])
  }
  #legend(leg.pos,colnames(x$dframe)[var.subset], col = catleg, pch = 22)
}

#----------------------------------end jointplot--------------------------------


}
