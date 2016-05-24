# R code from Remko Duursma (crownhull.R YplantQMC r package)
crownhull <- function(xyz, plotit=TRUE,col="forestgreen", alpha=0.8){
    
  # construct the hull (gives area and volume)
  if(is.list(xyz) && !is.data.frame(xyz))
    p <- as.matrix(do.call("rbind", xyz))
  else
    p <- as.matrix(xyz)
  
  ch <- convhulln(p, "FA")
  
  # construct the hull for plotting
  if(plotit){
    ch2 <- t(convhulln(p, "Qt"))
    triangles3d(p[ch2,1],p[ch2,2],p[ch2,3],col=col,alpha=alpha)
  }
  
  return(list(crownvolume=ch$vol, crownsurface=ch$area))
}