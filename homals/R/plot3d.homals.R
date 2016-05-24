`plot3d.homals` <-
function(x, plot.dim = c(1,2,3), plot.type = "jointplot", var.subset, type, xlab, ylab, zlab, col, main, sphere = TRUE, 
         bgpng = "particle.png", ax.grid = TRUE, ...)
{
#S3 plot method for objects of class "homals"
#Produces various 3D-plots
#plot.dim ... vector of length 3 with dimensions to be plotted against

if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !")
if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
pd1 <- plot.dim[1]
pd2 <- plot.dim[2]
pd3 <- plot.dim[3]
if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")


x1 <- x$objscores[,pd1]
y1 <- x$objscores[,pd2]
z1 <- x$objscores[,pd3]
nvar <- dim(x$dframe)[2]
if (missing(var.subset)) var.subset <- 1:nvar

if (missing(xlab)) xlab <- paste("Dimension",pd1)
if (missing(ylab)) ylab <- paste("Dimension",pd2)
if (missing(zlab)) zlab <- paste("Dimension",pd3)

if (is.null(bgpng)) {
    texture1 <- NULL
} else {
   texture1 <- system.file(paste("textures/",bgpng,sep=""), package = "rgl")
}

#------------------------------------objplot------------------------------------
#draws labeled object score plot
if (plot.type == "objplot") 
{
  
  if (missing(col)) col <- 4
  if (missing(main)) main1 <- "Object Plot"  else main1 = main
  
  #open3d()
  rgl.open()
  rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)
  text3d(x$objscores[,c(pd1,pd2,pd3)],texts = rownames(x$objscores), col = col, ...)
  axes3d(c('x','y','z'), labels = TRUE, color = "black", ...)
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black", ...)
  if (ax.grid) grid3d(c('x','y','z')) 
}
#---------------------------------end objplot-----------------------------------

#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- t(sapply(x$loadings, function(xy) xy[1,c(pd1,pd2,pd3)]))  #first solution only
  if (missing(main)) main1 <- "Loadings plot" else main1 <- main
 
  rgl.open()
  rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)
  points3d(xycoor, col = "RED", size = 2)
  text3d(xycoor, texts = rownames(xycoor), col = "RED")
  
  for (i in 1:nvar)
    lines3d(rbind(xycoor[i,],c(0,0,0)), col = 4)
    
  axes3d(c('x','y','z'), labels = TRUE, color = 1)
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black")
  if (ax.grid) grid3d(c('x','y','z'))
  
}
#-------------------------------- end loadplot ---------------------------------



#--------------------------------------labplot----------------------------------
#plot labeled object scores (for each variable separately)

if (plot.type == "labplot") {
  if (missing(col)) col <- 4
  for (i in var.subset) 
  {
    if (missing(main)) main1 <- paste("Label Plot for",colnames(x$dframe)[i])  else main1 = main
    
    rgl.open()
    rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)
    text3d(x1, y1, z1, texts = as.vector(x$dframe[,i]), col = col, ...)
    axes3d(c('x','y','z'), labels = TRUE, color = "black", ...)
    title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black", ...)
    if (ax.grid) grid3d(c('x','y','z'))
  }   
} 
#-----------------------------------end labplot---------------------------------

#----------------------------------catplot--------------------------------------
#plots the rank-restricted category quantifications for each variable

if (plot.type == "catplot") {

  if (missing(col)) col <- 4
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Category plot for",colnames(x$dframe[i]))  else main1 <- main
    
    rgl.open()
    rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)
    text3d(x$catscores[[i]][,c(pd1,pd2,pd3)], texts = levels(x$dframe[,i]), col = col, ...)
    axes3d(c('x','y','z'), labels = TRUE, color = "black", ...)
    title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black", ...)
    if (ax.grid) grid3d(c('x','y','z'))
  }
}
#----------------------------------end catplot----------------------------------

#---------------------------------- starplot------------------------------------

if (plot.type == "starplot") {
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Star plot for", colnames(x$dframe[i])) else main1 <- main
    y <- computeY(x$dframe[,i], x$objscores[,c(pd1,pd2,pd3)])
    rgl.open()
    rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)
    points3d(x1, y1, z1, col = 4, size = 3)
    text3d(y, texts = rownames(y), col = "RED")
    
    for (j in 1:length(x$dframe[,i]))                      #lineas
      lines3d(rbind(x$objscores[j,c(pd1,pd2,pd3)],y[x$dframe[,i][j],]), col = 4)
    
    axes3d(c('x','y','z'), labels = TRUE, color = "black")
    title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black")
    if (ax.grid) grid3d(c('x','y','z'))
  }     
}
#----------------------------------end starplot---------------------------------

#----------------------------------end jointplot--------------------------------

if (plot.type == "jointplot") {
   
  if (missing(main)) main1 <- "Joint Plot" else main1 <- main
  catcol <- rainbow(nvar)

  rgl.open()
  rgl.bg(sphere = sphere, texture = texture1, back = "filled", color = "white", ...)  
  text3d(x1, y1, z1, texts = rownames(x$dframe), col = 1) 
  for (i in var.subset) 
    text3d((x$catscores[[i]][,c(pd1,pd2,pd3)]), texts = rownames(x$catscores[[i]]), col = catcol[i])
  
  axes3d(c('x','y','z'), labels = TRUE, color = "darkgray")
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "darkgray")
  if (ax.grid) grid3d(c('x','y','z'))
}  
  
#----------------------------------end jointplot--------------------------------

 dev.off()



}

