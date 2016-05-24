####### R function: plot3D ###########

# For plotting points in 3D space.

# Last changed: 07 MAR 2016

plot3D <- function(x,y,z,main=NULL,
                   xlab=deparse(substitute(x)),
                   ylab=deparse(substitute(y)),
                   zlab=deparse(substitute(z)),
                   xlim=NULL,ylim=NULL,zlim=NULL,
                   cex=1,ptCol="red",axCol="navy",tickCol="navy",
                   labCol="DarkGreen",bgCol="white",addBase=FALSE,
                   baseCol="cyan",ptAlpha=0.7,baseAlpha=0.7,type="p")

{
   # Create working versions of input data:

   wx <- x ; wy <- y ; wz <- z

   # Remove data outside of specified ranges:

   omitInds <- NULL
   if (!is.null(xlim))
     omitInds <- union(omitInds,(1:length(x))[(wx<xlim[1])|(wx>xlim[2])])
   if (!is.null(ylim))
     omitInds <- union(omitInds,(1:length(y))[(wy<ylim[1])|(wy>ylim[2])])
   if (!is.null(zlim))
     omitInds <- union(omitInds,(1:length(z))[(wz<zlim[1])|(wz>zlim[2])])

   if (!is.null(omitInds))
   {
      wx <- wx[-omitInds] ; wy <- wy[-omitInds] ; wz <- wz[-omitInds]
   }
   
   # Insert safeguards against too much data being omitted:

   if (length(wx)==0) wx <- quantile(x,c(0.05,0.5,0.95))
   if (length(wy)==0) wy <- quantile(y,c(0.05,0.5,0.95))
   if (length(wz)==0) wz <- quantile(z,c(0.05,0.5,0.95))

   pwx <- wx ; pwy <- wy ; pwz <- wz

   if (!is.null(xlim)) pwx <- c(pwx,xlim) 
   if (!is.null(ylim)) pwy <- c(pwy,ylim)
   if (!is.null(zlim)) pwz <- c(pwz,zlim)

   xtcks <- pretty(pwx) ; ytcks <- pretty(pwy) ; ztcks <- pretty(pwz)

   xtcks <- xtcks[(xtcks>min(wx))&(xtcks<max(wx))]
   ytcks <- ytcks[(ytcks>min(wy))&(ytcks<max(wy))]
   ztcks <- ztcks[(ztcks>min(wz))&(ztcks<max(wz))]

   # Insert safeguards against too much data being omitted:
   
   if (length(xtcks)==0)
   {
     if (is.null(xlim)) 
        xtcks <- quantile(x,c(0.25,0.75))
     if (!is.null(xlim)) 
        xtcks <- xlim
   }
   if (length(ytcks)==0)
   {
     if (is.null(ylim)) 
        ytcks <- quantile(y,c(0.25,0.75))
     if (!is.null(ylim)) 
        ytcks <- ylim
   }
   if (length(ztcks)==0)
   {
     if (is.null(zlim)) 
        ztcks <- quantile(z,c(0.25,0.75))
     if (!is.null(zlim)) 
        ztcks <- zlim
   }
   
   # Set axis limits:
   
   a.x <- min(wx) ; b.x <- max(wx)   
   a.y <- min(wy) ; b.y <- max(wy)
   a.z <- min(wz) ; b.z <- max(wz)

   # Transform data to unit cube:
     
   tx <- tranUnitInt(wx,a.x,b.x)
   ty <- tranUnitInt(wy,a.y,b.y)
   tz <- tranUnitInt(wz,a.z,b.z)

   # Set up axes for unit cube:
   
   rgl.clear()
   rgl.bg(color=bgCol)
   rgl.lines3d(c(0,1.1),rep(0,2),rep(0,2),
                size=3,col=axCol,add=TRUE)
   rgl.lines3d(rep(0,2),c(0,1.1),rep(0,2),
               size=3,col=axCol,add=TRUE)
   rgl.lines3d(rep(0,2),rep(0,2),c(0,1.1),
                size=3,col=axCol,add=TRUE)

   rgl.texts(1.1,0,0,"x",col=axCol)
   rgl.texts(0,1.1,0,"y",col=axCol)
   rgl.texts(0,0,1.1,"z",col=axCol)
   
   rgl.texts(1.15,0,0,xlab,col=labCol)
   rgl.texts(0,1.15,0,ylab,col=labCol)
   rgl.texts(0,0,1.15,zlab,col=labCol)
   
   if (addBase)
      rgl.quads(c(-0.1,-0.1,1.1,1.1),
                 c(-0.1,1.1,1.1,-0.1),
                 rep(0,4),col=baseCol,alpha=baseAlpha)

  # Add `main' string:

  if (!is.null(main))
     rgl.texts(0.5,1.15,0,main,col=axCol)

  # Add tick marks:

  txtcks <- tranUnitInt(xtcks,a.x,b.x)
  for (ix in 1:length(txtcks))
  {
     rgl.lines3d(rep(txtcks[ix],2),rep(0,2),c(0,0.015),col=axCol,add=TRUE)
     rgl.lines3d(rep(txtcks[ix],2),c(0,0.015),rep(0,2),col=axCol,add=TRUE)
     rgl.texts(txtcks[ix],-0.05,-0.05,as.character(xtcks[ix]),col=tickCol)
  }
   
  tytcks <- tranUnitInt(ytcks,a.y,b.y)
  for (iy in 1:length(tytcks))
  {
     rgl.lines3d(rep(0,2),rep(tytcks[iy],2),c(0,0.015),col=axCol,add=TRUE)
     rgl.lines3d(c(0,0.015),rep(tytcks[iy],2),rep(0,2),col=axCol,add=TRUE)
     rgl.texts(-0.05,tytcks[iy],-0.05,as.character(ytcks[iy]),col=tickCol)
  }

  tztcks <- tranUnitInt(ztcks,a.z,b.z)
  for (iz in 1:length(tztcks))
  {
     rgl.lines3d(rep(0,2),c(0,0.015),rep(tztcks[iz],2),col=axCol,add=TRUE)
     rgl.lines3d(c(0,0.015),rep(0,2),rep(tztcks[iz],2),col=axCol,add=TRUE)
     rgl.texts(-0.05,-0.05,tztcks[iz],as.character(ztcks[iz]),col=tickCol)
  }

  # Add points:
   
  if (type=="p")
     rgl.spheres(tx,ty,tz,col=ptCol,radius=cex*0.015)
################
#  CHANGE FOR 12 AUG 2010 TALK   
#     rgl.spheres(tx,ty,tz,col=ptCol,alpha=ptAlpha,radius=cex*0.015)
###########################
   
  # Add limit information as attributes:
   
  obj <- 0
  attr(obj,"a.x") <- a.x ; attr(obj,"b.x") <- b.x
  attr(obj,"a.y") <- a.y ; attr(obj,"b.y") <- b.y
  attr(obj,"a.z") <- a.z ; attr(obj,"b.z") <- b.z

  invisible(obj)
}

############# End of plot3D ##############

