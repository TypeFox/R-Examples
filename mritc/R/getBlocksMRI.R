getBlocksMRI <- function(mask, nblocks){
    if(length(dim(mask)) != 3)
        stop("The 'mask' has to be of dimension 3.")
    if(! all(unique(as.vector(mask)) %in% 0:1))
        stop("The value of mask has to be either 1 or 0.")
    if(all(mask==0))
        stop("All voxels are outside the mask.")
    if(! nblocks %in% c(2,8))
        stop("The number of blocks has to be equal to either 2 or 8.")
    
    
    ind <- which(mask==1, arr.ind=T)
    x.mark <- ind[,1] %% 2
    y.mark <- ind[,2] %% 2
    z.mark <- ind[,3] %% 2
    
    if(nblocks==2){
        BW <- ifelse(rowSums(cbind(x.mark, y.mark, z.mark)) %% 2 == 1, 1, 0)
        lapply(1:2, function(x) which(BW == (x-1)))
    }
    
    else
      lapply(1:8, function(x) which((x.mark + y.mark*2 + z.mark*4) == (x-1)))
}

