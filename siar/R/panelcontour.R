panelcontour <-
function(x, y, ...)
{
     usr <- par("usr"); on.exit(par(usr))
     par(usr = c(usr[1:2], 0, 1.5) )
     kd <- kde2d(x,y)
     kdmax <- max(kd$z)
     #points(x,y,col="light grey")
     contour(kd,add=TRUE,drawlabels=FALSE,levels=c(kdmax*0.01,kdmax*0.5,kdmax*0.9),col=c("red","blue","green"),...)
}
