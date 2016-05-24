assign("network.design",
function(formula, vgm.model, xmin, xmax, ymin, ymax, npoint.x, npoint.y, npoints, boundary=NULL, type, ...){
if (is.null(boundary)) {
grid.p<-expand.grid(x=seq(xmin,xmax,(xmax-xmin)/(npoint.x-1)), y=seq(ymin,ymax,(ymax-ymin)/(npoint.y-1)))
plot(grid.p,pch=19,cex=0.5)
grid.p$z <- grid.p$x
}
else if (is.null(boundary)==FALSE) {
df.pts<-spsample(boundary, n=npoints, type="regular")
plot(boundary)
points(df.pts,pch=19,cex=0.5)
grid.p <- data.frame(df.pts,df.pts@coords[,1])
names(grid.p) <- c("x","y","z")
}
K = krige.cv(formula=formula, ~x+y, grid.p, vgm.model, ...)
ASEPE <- mean((K[,2])^0.5)              
ASEPE
}
)
