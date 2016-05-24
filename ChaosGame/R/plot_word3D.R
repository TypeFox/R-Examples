plot_word3D <- 
function(word = "copula", R = 20, plot.rgl = TRUE, copula = TRUE, portion = 0.2, color.rgl.plot = "green2magenta", plot.surface = "Sphere", histogram = TRUE, shift = 1.2, orbit = 3000, cex.label = 0.7, size.lines = 0.1){
  
  A <- .chaos_game_word(word = word, R = R, shift = shift, orbit = orbit)
  
  if(plot.surface == "Sphere"){surface <- .plot_ball(A)}
  if(plot.surface == "Helix"){surface <-  .plot_Helix(A)}
  if(plot.surface == "Torus"){surface <- .plot_Torus(A)}
  if(plot.surface == "EnneperMinimalSurface"){surface <- .plot_EnneperMinimalSurface(A)}
  if(plot.surface == "CatalanSurface"){surface <- .plot_CatalanSurface(A)}
  
  if(copula == TRUE){
    
    perc<-round(nrow(surface)*portion)
    trafox<-ecdf(surface$x[sample(seq(1:length(surface$z)), perc, replace = FALSE)])
    trafoy<-ecdf(surface$y[sample(seq(1:length(surface$z)), perc, replace = FALSE)])
    trafoz<-ecdf(surface$z[sample(seq(1:length(surface$z)), perc, replace = FALSE)])
    surface<-data.frame(x=trafox(surface$x),y=trafoy(surface$y),z=trafoz(surface$z))
    
  }else{
    
    x <- (surface$x-min(surface$x))/(max(surface$x)-min(surface$x))
    y <- (surface$y-min(surface$y))/(max(surface$y)-min(surface$y))
    z <- (surface$z-min(surface$z))/(max(surface$z)-min(surface$z))
    surface<-data.frame(x=x,y=y,z=z)
    
  }
  
  hx <- hist(surface[,2], plot=FALSE)
  hxs <- hx$density 
  hy <- hist(surface[,1], plot=FALSE)
  hys <- hy$density 
  hz <- hist(surface[,3], plot=FALSE)
  hzs <- hz$density
  
  hxs_max <- 1
  hys_max <- 1 
  hzs_max <- 1 
  
  if(copula == FALSE){
    
    hxs_max <- max(hx$density)
    hys_max <- max(hy$density) 
    hzs_max <- max(hz$density) 
    
    hxs <- hx$density/hxs_max
    hys <- hy$density/hys_max
    hzs <- hz$density/hzs_max
    
  }
  
  ## no overlap in the adjoining corner
  xmax <- tail(hx$breaks, n=1) + diff(tail(hx$breaks, n=2))
  ymax <- tail(hy$breaks, n=1) + diff(tail(hy$breaks, n=2))
  zmax <- tail(hz$breaks, n=1) + diff(tail(hz$breaks, n=2))
  
  dd <- 0.75
  col_lines <- "gray5"

  
  
  
  # rgl plot or scatter3D
  
  if(plot.rgl){
    
    col_hist <- "gray90"
    
    od<-order(surface$z)
    df_sph_coord<-car2sph(surface$x[od], surface$y[od], surface$z[od], deg = FALSE)
    deg = 15
    open3d() 
    rgl.sphpoints(df_sph_coord, deg=FALSE, col=.colorfunction(length(surface$z), col = color.rgl.plot), size=0.1)
    
    if(histogram){
      
      # Box Histogramm x y and z:
      plot3d(data.frame(x=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2,y=c(-0.15,-0.15,-0.65,-0.65,-0.15)*1.1+0.03,z=c(-0.05,-0.05,-0.05,-0.05,-0.05)),col="gray50",add=T,type='l') 
      plot3d(data.frame(x=c(1.05,1.05,1.55,1.55,1.05)*1.1-0.02,y=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2,z=c(-0.05,-0.05,-0.05,-0.05,-0.05)),col="gray50",add=T,type='l') 
      plot3d(data.frame(x=c(-0.05,-0.05,-0.05,-0.05,-0.05),y=c(-0.15,-0.15,-0.65,-0.65,-0.15)*1.1+0.03,z=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2),col="gray50",add=T,type='l') 
      text3d(c(1.05*dd+(1-dd)/2,1.05*dd+(1-dd)/2,1.15,1.65,1.15,1.65,-0.05,-0.05),c(-0.15,-0.65,-0.05*dd+(1-dd)/2,-0.05*dd+(1-dd)/2,1.07*dd+(1-dd)/2,1.07*dd+(1-dd)/2,-0.15,-0.65),c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05,1.07*dd+(1-dd)/2,1.07*dd+(1-dd)/2),c("0",paste(round(hxs_max,2)),"0",paste(round(hys_max,2)),"0",paste(round(hys_max,2)),"0",paste(round(hzs_max,2))),cex=0.5,col="gray50") 
      text3d(c(-0.07*dd+(1-dd)/2,-0.07*dd+(1-dd)/2,-0.05,-0.05),c(-0.15,-0.65,-0.15,-0.65),c(-0.05,-0.05,-0.05*dd+(1-dd)/2,-0.05*dd+(1-dd)/2),c("0",paste(round(hxs_max,2)),"0",paste(round(hzs_max,2))),cex=0.5,col="gray50") 
      
      ## manually create each histogram
      for (ii in seq_along(hx$counts)) {
        x <- (hx$breaks[ii]*c(.9,.9,.1,.1) + hx$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
        y <- -hxs[ii]*c(0,0.5,0.5,0)-0.15
        z <- rep(ymax, 4)-1.1
        quads3d(x,y,z, color=col_hist, alpha=.6)
        plot3d(data.frame(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1])),col="gray20",add=T,type='l') 
      }
      
      for (ii in seq_along(hy$counts)) {
        x <- hys[ii]*c(0,0.5,0.5,0)+1.15
        y <- (hy$breaks[ii]*c(.9,.9,.1,.1) + hy$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
        z <- rep(xmax, 4)-1.1
        quads3d(x,y,z, color=col_hist, alpha=.6)
        plot3d(data.frame(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1])),col="gray20",add=T,type='l') 
      }
      
      for (ii in seq_along(hz$counts)) {
        x <- rep(ymax, 4)-1.1
        y <- -hzs[ii]*c(0,0.5,0.5,0)-0.15
        z <- (hz$breaks[ii]*c(.9,.9,.1,.1) + hz$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
        quads3d(x,y,z, color=col_hist, alpha=.6)
        plot3d(data.frame(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1])),col="gray20",add=T,type='l') 
      }
    }
    
    
    # Box:
    plot3d(data.frame(x=c(0,1,1,0,0,0,1,1,0,0)*1.1-0.05,y=c(0,0,1,1,0,0,0,1,1,0)*1.1-0.05,z=c(0,0,0,0,0,1,1,1,1,1)*1.1-0.05),col=col_lines,add=T,type='l') 
    plot3d(data.frame(x=c(1,1)*1.1-0.05,y=c(0,0)*1.1-0.05,z=c(0,1)*1.1-0.05),col=col_lines,add=T,type='l') 
    plot3d(data.frame(x=c(0,0)*1.1-0.05,y=c(1,1)*1.1-0.05,z=c(0,1)*1.1-0.05),col=col_lines,add=T,type='l') 
    plot3d(data.frame(x=c(1,1)*1.1-0.05,y=c(1,1)*1.1-0.05,z=c(0,1)*1.1-0.05),col=col_lines,add=T,type='l') 
    
    # Label:
    text3d(c(0.5,1,0),c(0,0.5,0),c(0,0,0.5),c("x","y","z"),cex=1,col=col_lines)
    text3d(c(0,0.5,1)*dd+(1-dd)/2,c(1.05,1.05,1.05),c(1.12,1.12,1.12),c("0","0.5","1"),cex=0.7,col=col_lines) # x Achse
    text3d(c(-0.05,-0.05,-0.05),c(0,0.5,1)*dd+(1-dd)/2,c(1.12,1.12,1.12),c("0","0.5","1"),cex=0.7,col=col_lines) # y Achse
    text3d(c(1.12,1.12,1.12),c(1.12,1.12,1.12),c(0,0.5,1)*dd+(1-dd)/2,c("0","0.5","1"),cex=0.7,col=col_lines) # z Achse
    plot3d(data.frame(x=c(1.05,1.07),y=c(1.05,1.07),z=c(0,0)*dd+(1-dd)/2),col=col_lines,add=T,type="l")  # z Achse
    plot3d(data.frame(x=c(1.05,1.07),y=c(1.05,1.07),z=c(0.5,0.5)*dd+(1-dd)/2),col=col_lines,add=T,type="l")  # z Achse
    plot3d(data.frame(x=c(1.05,1.07),y=c(1.05,1.07),z=c(1,1)*dd+(1-dd)/2),col=col_lines,add=T,type="l")  # z Achse
    
    plot3d(data.frame(x=c(0.5,0.5)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # x Achse
    plot3d(data.frame(x=c(0,0)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # x Achse
    plot3d(data.frame(x=c(1,1)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # x Achse
    
    plot3d(data.frame(x=c(-0.05,-0.05),y=c(0.5,0.5)*dd+(1-dd)/2,z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # y Achse
    plot3d(data.frame(x=c(-0.05,-0.05),y=c(0,0)*dd+(1-dd)/2,z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # y Achse
    plot3d(data.frame(x=c(-0.05,-0.05),y=c(1,1)*dd+(1-dd)/2,z=c(1.05,1.07)),col=col_lines,add=T,type="l")  # y Achse
    
    plot3d(data.frame(x=surface$x*dd+(1-dd)/2,y=surface$y*dd+(1-dd)/2,z=rep(-0.05,length(surface$y))),col='grey85',add=T,type='p', size = 0.001)
    plot3d(data.frame(x=surface$x*dd+(1-dd)/2,y=rep(1.05,length(surface$y)),z=surface$z*dd+(1-dd)/2),col='grey85',add=T,type='p', size = 0.001)
    plot3d(data.frame(x=rep(-0.05,length(surface$y)),y=surface$y*dd+(1-dd)/2,z=surface$z*dd+(1-dd)/2),col='grey85',add=T,type='p', size = 0.001)
    
    
  }else{
    
    color <- ramp.col(c("chartreuse","green","darkolivegreen","magenta","deeppink4","darkmagenta")) 
    #color <- ramp.col(c("darkmagenta","deeppink4","magenta","darkolivegreen","green","chartreuse")) 
    col_hist <- "darkslategray3" 
    col_lines <- "gray5"
    
    # Box:
    scatter3D(x=c(0,1,1,0,0,0,1,1,0,0)*1.1-0.05,y=c(0,0,1,1,0,0,0,1,1,0)*1.1-0.05,z=c(0,0,0,0,0,1,1,1,1,1)*1.1-0.05,col=col_lines,type='l', xlim=c(0,1.3), ylim=c(-0.3,1), zlim=c(-0.3,1.1),
              pch = 19, cex = 0.1, xlab = "", ylab = "", 
              zlab = "", clab = c("",""),
              main = "", ticktype = "detailed", box = FALSE,
              theta = 40, d = 2, phi = 30,   #  theta gives the azimuthal direction and phi the colatitude
              colkey = FALSE,lwd=size.lines) 
    
    scatter3D(x=c(1,1)*1.1-0.05,y=c(0,0)*1.1-0.05,z=c(0,1)*1.1-0.05,col=col_lines,add=T,type='l',lwd=size.lines) 
    scatter3D(x=c(0,0)*1.1-0.05,y=c(1,1)*1.1-0.05,z=c(0,1)*1.1-0.05,col=col_lines,add=T,type='l',lwd=size.lines) 
    scatter3D(x=c(1,1)*1.1-0.05,y=c(1,1)*1.1-0.05,z=c(0,1)*1.1-0.05,col=col_lines,add=T,type='l',lwd=size.lines) 
    
    # Box Histogramm x y und z:
    scatter3D(x=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2,y=c(-0.15,-0.15,-0.65,-0.65,-0.15)*1.1+0.03,z=c(-0.05,-0.05,-0.05,-0.05,-0.05),col="gray10",add=T,type='l',lwd=size.lines) 
    scatter3D(x=c(1.05,1.05,1.55,1.55,1.05)*1.1-0.02,y=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2,z=c(-0.05,-0.05,-0.05,-0.05,-0.05),col="gray10",add=T,type='l',lwd=size.lines) 
    scatter3D(x=c(-0.05,-0.05,-0.05,-0.05,-0.05),y=c(-0.15,-0.15,-0.65,-0.65,-0.15)*1.1+0.03,z=(c(0,1,1,0,0)*1.05-0.025)*dd+(1-dd)/2,col="gray10",add=T,type='l',lwd=size.lines) 
    
    text3D(c(1.05*dd+(1-dd)/2,1.05*dd+(1-dd)/2,1.15,1.65,1.15,1.65,-0.05,-0.05),c(-0.15,-0.65,-0.05*dd+(1-dd)/2,-0.05*dd+(1-dd)/2,1.07*dd+(1-dd)/2,1.07*dd+(1-dd)/2,-0.15,-0.65),c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05,1.07*dd+(1-dd)/2,1.07*dd+(1-dd)/2),c("0",paste(round(hxs_max,2)),"0",paste(round(hys_max,2)),"0",paste(round(hys_max,2)),"0",paste(round(hzs_max,2))),cex=cex.label*0.7,col="gray10",add=T) 
    text3D(c(-0.07*dd+(1-dd)/2,-0.07*dd+(1-dd)/2,-0.05,-0.05),c(-0.15,-0.65,-0.15,-0.65),c(-0.05,-0.05,-0.05*dd+(1-dd)/2,-0.05*dd+(1-dd)/2),c("0",paste(round(hxs_max,2)),"0",paste(round(hzs_max,2))),cex=cex.label*0.7,col="gray10",add=T) 
    
    
    # Label:
    text3D(c(0.5,1,0),c(0,0.5,0),c(0,0,0.5),c("x","y","z"),cex=cex.label*1.4,col=col_lines,add=T)
    text3D(c(0,0.5,1)*dd+(1-dd)/2,c(1.05,1.05,1.05),c(1.12,1.12,1.12),c("0","0.5","1"),cex=cex.label,col=col_lines,add=T) # x Achse
    text3D(c(-0.05,-0.05,-0.05),c(0,0.5,1)*dd+(1-dd)/2,c(1.12,1.12,1.12),c("0","0.5","1"),cex=cex.label,col=col_lines,add=T) # y Achse
    text3D(c(1.12,1.12,1.12),c(1.12,1.12,1.12),c(0,0.5,1)*dd+(1-dd)/2,c("0","0.5","1"),cex=cex.label,col=col_lines,add=T) # z Achse
    scatter3D(x=c(1.05,1.07),y=c(1.05,1.07),z=c(0,0)*dd+(1-dd)/2,col=col_lines,add=T,type="l",lwd=size.lines)  # z Achse
    scatter3D(x=c(1.05,1.07),y=c(1.05,1.07),z=c(0.5,0.5)*dd+(1-dd)/2,col=col_lines,add=T,type="l",lwd=size.lines)  # z Achse
    scatter3D(x=c(1.05,1.07),y=c(1.05,1.07),z=c(1,1)*dd+(1-dd)/2,col=col_lines,add=T,type="l",lwd=size.lines)  # z Achse
    
    scatter3D(x=c(0.5,0.5)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # x Achse
    scatter3D(x=c(0,0)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # x Achse
    scatter3D(x=c(1,1)*dd+(1-dd)/2,y=c(1.05,1.05),z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # x Achse
    
    scatter3D(x=c(-0.05,-0.05),y=c(0.5,0.5)*dd+(1-dd)/2,z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # y Achse
    scatter3D(x=c(-0.05,-0.05),y=c(0,0)*dd+(1-dd)/2,z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # y Achse
    scatter3D(x=c(-0.05,-0.05),y=c(1,1)*dd+(1-dd)/2,z=c(1.05,1.07),col=col_lines,add=T,type="l",lwd=size.lines)  # y Achse
    
    
    ## manually create each histogram
    for (ii in seq_along(hx$counts)) {
      x <- (hx$breaks[ii]*c(.9,.9,.1,.1) + hx$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
      y <- -hxs[ii]*c(0,0.5,0.5,0)-0.15
      z <- rep(ymax, 4)-1.1
      polygon3D(x,y,z, color=col_hist, alpha=.6,add=T)
      scatter3D(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1]),col="gray20",add=T,type='l',lwd=size.lines) 
    }
    
    
    for (ii in seq_along(hy$counts)) {
      x <- hys[ii]*c(0,0.5,0.5,0)+1.15
      y <- (hy$breaks[ii]*c(.9,.9,.1,.1) + hy$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
      z <- rep(xmax, 4)-1.1
      polygon3D(x,y,z, color=col_hist, alpha=.6,add=T)
      scatter3D(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1]),col="gray20",add=T,type='l',lwd=size.lines) 
    }
    
    for (ii in seq_along(hz$counts)) {
      x <- rep(ymax, 4)-1.1
      y <- -hzs[ii]*c(0,0.5,0.5,0)-0.15
      z <- (hz$breaks[ii]*c(.9,.9,.1,.1) + hz$breaks[ii+1]*c(.1,.1,.9,.9))*dd+(1-dd)/2
      polygon3D(x,y,z, color=col_hist, alpha=.6,add=T)
      scatter3D(x=c(x,x[1]),y=c(y,y[1]),z=c(z,z[1]),col="gray20",add=T,type='l',lwd=size.lines) 
    }
    
    scatter3D(x=c(rep(-0.05,length(surface$y)),surface$x*dd+(1-dd)/2,surface$x*dd+(1-dd)/2), y=c(surface$y*dd+(1-dd)/2,rep(1.05,length(surface$y)),surface$y*dd+(1-dd)/2), z=c(surface$z*dd+(1-dd)/2,surface$z*dd+(1-dd)/2,rep(-0.05,length(surface$y))), add = TRUE, col = 'grey85',
              pch = 19, cex = 0.1, ticktype = "detailed")
    
    scatter3D(x = surface$x*dd+(1-dd)/2, y = surface$y*dd+(1-dd)/2, z = surface$z*dd+(1-dd)/2, xlim=c(0,1.3), ylim=c(-0.3,1), zlim=c(-0.3,1), col = color,
              pch = 19, cex = 0.1, xlab = "", ylab = "", add = TRUE,
              zlab = "", clab = c("",""),
              main = "", ticktype = "detailed", box = FALSE,
              theta = 40, d = 2, phi = 30,   
              colkey = FALSE)
    
  }
  
  return(surface)
  
}