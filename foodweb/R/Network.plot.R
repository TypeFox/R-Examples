
create.plot <- function(foodweb, radii, cols, levels) {

library(rgl)

if (missing(cols)==TRUE) { 
  cols <- as.vector(c("red", "blue", "green", "purple", "black", "orange", "pink", "grey", "brown", "turquoise"))
  cols <- cols[1:levels]}

# Transform 360 degrees into radians
deg <- 360*pi/180
S <- ncol(foodweb)

#Create a matrix to store coordinates of the points and their colour
net.plot <- matrix(nrow= S, ncol=4)
colnames(net.plot) <- c("x", "y","z", "colour")
net.plot[,2] <- t(foodweb[S+1,])
colour.code <- cbind(cols)
rownames(colour.code) <- unique(t((foodweb[S+1,])))

# For setting radius at each trophic level
if (exists("radii")) {
  radius.code <- matrix(ncol=length(unique(t((foodweb[S+1,])))))
  radius.code <- cbind(radii)
  rownames(radius.code) <- unique(t((foodweb[S+1,])))
  } else {
  for (level in unique(t((foodweb[S+1,])))) {
      n <- length(which(net.plot[,2]==level))
      radius.code[which(colnames(radius.code)==level)] <- n
    }
  }

#For calculation of the coordinates of each species
for (level in unique(t((foodweb[S+1,])))) {
  n <- length(which(net.plot[,2]==level))
  r <- as.double(radius.code[which(rownames(radius.code)==level)])
  x <-as.vector(r*(cos(seq(0, deg, len = n+1))))
  z <-as.vector(r*(sin(seq(0, deg, len = n+1))))
  
  net.plot[which(net.plot[,2]==level),1] <- x[-(length(x))]
  net.plot[which(net.plot[,2]==level),3] <- z[-(length(z))]
  net.plot[which(net.plot[,2]==level),4] <- as.character(colour.code[which(rownames(colour.code)==level)])
}

#Place spheres in space
plot3d(net.plot[,1],net.plot[,3], net.plot[,2], type="s", col=net.plot[,4], size=2, box=FALSE, axes=FALSE, xlab="", ylab="", zlab="", aspect=TRUE, top=TRUE)

#Add lines where trophic links exist
for (sp in 1:S) {
  for (prey in intersect(which(foodweb[1:S,sp]==1), which(foodweb[S+1,]!=foodweb[S+1,sp]))) {
    lines3d(c(net.plot[sp,1],net.plot[prey,1]), c(net.plot[sp,3],net.plot[prey,3]), c(net.plot[sp,2],net.plot[prey,2]), cex=0.1)
  }
}


}