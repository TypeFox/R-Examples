##
## onle for debugging purposes: visualize RBF surrogate (only d=2)
## 
## caller: cobraPhaseII.R, switch DEBUG_RBF 
##
drawSurrogate3d <- function(surrogate,fName,X,Y,lower,upper,add=FALSE,newWindow=TRUE) {
  if (!is.null(surrogate)) {
    fac <- (upper[1]-lower[1])/2
    N = nrow(X)
    ZS = sapply(1:N,function(j)predict(surrogate,cbind(X[,j]/fac,Y[,j]/fac)))
    jet.colors <-   
      grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorzjet <- jet.colors(100)  # 100 separate colors 
    
    graphics::filled.contour(X[1,],Y[,1],ZS,main=sprintf("%s, surrogate, N=%d",fName,surrogate$npts),
                   color.palette=jet.colors)
    if (newWindow) {
      rgl::open3d()      
    } else {
      if (!add) rgl::rgl.clear()
    }
    rngZS=max(ZS)-min(ZS)
    Z2=ZS
    rngX = upper[1] - lower[1]
    rngY = upper[2] - lower[2]
    rngZ = max(ZS)-min(ZS)
    rgl::aspect3d(1/rngX,1/rngY,1/rngZ)
    rgl::rgl.surface(x=X, y=Y, z=Z2,
                     coords=c(1,2,3), 
                     color=colorzjet[ findInterval(Z2, seq(min(Z2), max(Z2), length=100))]
    )
    rgl::axes3d(color="black") # draw axes and box with tickmarks
    rgl::title3d(sprintf("%s, surrogate, N=%d",fName,surrogate$npts), '', 'x', 'y', 'z',color="black")
    rgl::rgl.bg(color="grey")
    rgl::rgl.bringtotop()
    return(ZS)
  }  
}
