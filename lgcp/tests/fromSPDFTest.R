library(lgcp)
library(spatstat)
library(sp)

polylist <- list()
# a pentagon
vert1 <- t(sapply((0:4)*(2*pi/5),function(theta){return(5*c(cos(theta),sin(theta)))}))
vert1 <- rbind(vert1,vert1[1,])
polylist[[1]] <- Polygons(list(Polygon(vert1)),ID="1")

# a diamond shape
vert2 <- matrix(c(-8,-10,-12,-10,0,2,0,-2),4,2)
vert2 <- rbind(vert2,vert2[1,])
polylist[[2]] <- Polygons(list(Polygon(vert2)),ID="2")

# a square with a rectangular hole
vert3 <- matrix(c(0,0,2,2,-10,-12,-12,-10),4,2)
vert3 <- rbind(vert3,vert3[1,])
vert4 <- matrix(c(0.25,0.75,0.75,0.25,-10.5,-10.5,-11.5,-11.5),4,2)
vert4 <- rbind(vert4,vert4[1,])
polylist[[3]] <- Polygons(list(Polygon(vert3),Polygon(vert4,hole=TRUE)),ID="3") 


sP <- SpatialPolygons(polylist)

df <- data.frame(list(atrisk=c(1,2,3)))

spdf <- SpatialPolygonsDataFrame(sP,df)

sar <- spatialAtRisk(spdf)
sar

rotsar <- affine(sar,mat=rotmat(pi/6))

# Now test simulation code:
#sim <- lgcpSim(affine(owin(poly=list(x=c(-15,-15,5,5),y=c(5,-15,-15,5))),rotmat(pi/6)),model.parameters=lgcppars(sigma=2,phi=1,theta=3),cellwidth=0.2,spatial.intensity=affine(sar,rotmat(pi/6)),temporal.intensity=function(x){100},plot=TRUE)
#sim <- lgcpSim(owin(poly=list(x=c(-15,-15,5,5),y=c(5,-15,-15,5))),model.parameters=lgcppars(sigma=2,phi=1,theta=3),cellwidth=0.2,spatial.intensity=sar,temporal.intensity=function(x){100},plot=TRUE)

