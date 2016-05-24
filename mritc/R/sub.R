#split on voxel into 8 subvoxels
#input is a mask
split18 <- function(mask){
  M <- array(rep(1,8), dim=rep(2,3))
  kronecker(mask, M)

}

#get the match up table of voxels and their subvoxels
#input
#  x: the mask
getPos18 <- function(mask){
  submask <- split18(mask==1)
  subpos <- which(submask==1)
  ind <- which(mask==1, arr.ind=T)
  d <- dim(mask)
  nx=d[1]; ny=d[2]; nz=d[3]
  v1 <- (ind[,3]-1)*8*nx*ny+(ind[,2]-1)*4*nx+(ind[,1]-1)*2+1
  v5 <- v1 + nx*ny*4
  pos <- cbind(v1, v1+1, v1+nx*2, v1+nx*2+1, v5, v5+1, v5+nx*2, v5+nx*2+1)
  matrix(match(as.vector(pos), subpos),  ncol=8)
}

