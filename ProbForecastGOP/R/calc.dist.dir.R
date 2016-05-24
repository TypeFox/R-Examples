"calc.dist.dir" <-
function(coord1,coord2,id,tol.angle.rad1,tol.angle.rad2,type){
  radius <- 6378.1
  mat.id <- matrix(cbind(coord1,coord2,id),ncol=3)
  unique.mat.id <- unique(mat.id)
  coord1.u <- as.numeric(unique.mat.id[,1])
  coord2.u <- as.numeric(unique.mat.id[,2])
  n.id <- length(coord1.u)
  dist.matrix <- matrix(0,nrow=n.id,ncol=n.id)  
  coord1.norm <- coord1.u/57.29577951
  coord2.norm <- coord2.u/57.29577951
  coord1.mat1 <- 
matrix(rep(coord1.norm,each=n.id),nrow=n.id,ncol=n.id,byrow=TRUE)
  coord1.mat2 <- t(coord1.mat1)
  coord2.mat1 <- 
matrix(rep(coord2.norm,each=n.id),nrow=n.id,ncol=n.id,byrow=TRUE)
  coord2.mat2 <- t(coord2.mat1)
  diff.coord1 <- coord1.mat1-coord1.mat2
  index.diff <- seq(1:(n.id^2))[diff.coord1 > pi]
  index.zero <- seq(1:(n.id^2))[diff.coord1==0]
  r.diff.coord1 <- replace(diff.coord1,index.diff,(2*pi)-diff.coord1[index.diff])
  sin.c2m1 <- sin(coord2.mat1)
  sin.c2m2 <- sin(coord2.mat2)
  cos.c2m1 <- cos(coord2.mat1)
  cos.c2m2 <- cos(coord2.mat2)
  cos.diff <- cos(r.diff.coord1)
  cos.A <- sin.c2m1*sin.c2m2+cos.c2m1*cos.c2m2*cos.diff
  cos.A2 <- cos.A*cos.A
  index.one <- seq(1:(n.id^2))[cos.A2 > 1]
  r.cos.A2 <- replace(cos.A2,index.one,1)
  sin.A <- sqrt(1-r.cos.A2)
  cos.angle.mat <- (sin.c2m2-(sin.c2m1*cos.A))/(cos.c2m1*sin.A)
  index.1 <- seq(1:(n.id^2))[cos.angle.mat > 1]
  index.m1 <- seq(1:(n.id^2))[cos.angle.mat < -1]
  r.cos.angle.mat <- replace(cos.angle.mat,index.1,1)
  r.cos.angle.mat <- replace(r.cos.angle.mat,index.m1,-1)
  angle.mat <- acos(r.cos.angle.mat)  # this is to determine the angle between the locations
                                    # in a matrix form
  r.angle.mat <- replace(angle.mat,index.zero,pi)
  if(type=="E"){
   index.angle <- seq(1:(n.id^2))[angle.mat >= tol.angle.rad1 & angle.mat <= tol.angle.rad2]
  }
  if(type=="N"){
   index.angle <- seq(1:(n.id^2))[angle.mat < tol.angle.rad1 | angle.mat > tol.angle.rad2]
  }
  index <- seq(1:(n.id^2))[cos.A>1]
  repl.cos.A <- replace(cos.A,index,1)
  dist.matrix <- round(radius*acos(repl.cos.A),2)
  dist.matrix.dir <- replace(dist.matrix,index.angle,0)
  
  return(dist.matrix.dir)
}
