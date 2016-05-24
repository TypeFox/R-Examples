"calc.dist" <-
function(coord1,coord2,id){
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
  sin.c2m1 <- sin(coord2.mat1)
  sin.c2m2 <- sin(coord2.mat2)
  cos.c2m1 <- cos(coord2.mat1)
  cos.c2m2 <- cos(coord2.mat2)
  cos.diff <- cos(diff.coord1)
  cos.A <- sin.c2m1*sin.c2m2+cos.c2m1*cos.c2m2*cos.diff
  index <- seq(1:(n.id^2))[cos.A>1]
  repl.cos.A <- replace(cos.A,index,1)
  dist.matrix <- round(radius*acos(repl.cos.A),2)
  
  return(dist.matrix)
}
