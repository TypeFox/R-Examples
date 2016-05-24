#split one voxel into 8 subvoxels
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


makeMRIspatial <- function(mask, nnei, sub=FALSE, bias=FALSE){
    if(! nnei %in% c(6, 18, 26))
        stop("'nnei' has to be among (6, 18, 26)")
    if(nnei %in% c(4, 6)) nblocks <- 2
    else nblocks <- 8

    weights <- NULL
    weineighbors <- NULL

    if(sub==FALSE){
        neighbors <- getNeighborsMRI(mask, nnei)
        blocks <- getBlocksMRI(mask, nblocks)
        subvox=NULL
 
    }
    else{
        submask <- split18(mask)
        neighbors <- getNeighborsMRI(submask, nnei)
        blocks <- getBlocksMRI(submask, nblocks)
        subvox <- getPos18(mask)
    }
    
    if(bias==TRUE){
        weights <-getWeightsMRI(nnei=nnei, sigma=1)
        nvertex <- nrow(neighbors)
        nei <- (neighbors < nvertex+1)
        nei <- nei %*% weights
        weineighbors <-  rowSums(nei) 
    }
    
    result <- list(neighbors=neighbors, blocks=blocks, sub=sub, subvox=subvox,
                   weights=weights, weineighbors=weineighbors)
 
    result
}

