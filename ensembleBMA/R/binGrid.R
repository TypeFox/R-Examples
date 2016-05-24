`binGrid` <-
function( val, coord1, coord2, nGrid = 65){
#
# copyright 2006-present, University of Washington. All rights reserved. 
# for terms of use, see the LICENSE file
#
   size <- nGrid
   grid <- matrix(NA, size,size)  
   paco <- cbind(coord1, coord2, val)
   x <- seq(min(coord1,na.rm=TRUE),max(coord1,na.rm=TRUE),length = size + 
1)  
   y <- seq(min(coord2,na.rm=TRUE),max(coord2,na.rm=TRUE),length = size + 
1)  
   k <- 1  
   while(k <= size){
      temp <- matrix(paco[paco[,2] >= y[k],], ncol = 3)
      temp <- matrix(temp[temp[,2] <= y[k+1],], ncol = 3)
      r <- 1
         while( r <= size){
           temmp <- matrix(temp[temp[,1] >= x[r],], ncol = 3)
           temmp <- matrix(temmp[temmp[,1] <= x[r+1],], ncol = 3)
        
           if(length(temmp) == 3){
              grid[size + 1 - k, r] <- temmp[,3]}
           if(length(temmp) > 3){
              grid[size + 1 - k, r] <- mean(temmp[,3])}
         r <- r + 1}
   k <- k + 1}
   
   temp <- grid
   k <- 1
   while(k <= size){  
      r <- 1
      while( r <= size){
        grid[r,size + 1 - k] <- temp[k,r]
        r <- r + 1}
   k <- k + 1}
  index.NA <- NULL
  for(i in 1:size){
    for(j in 1:size){ 
      if(is.na(grid[i,j])==TRUE){
        index.NA <- rbind(index.NA,c(i,j))
      }
    }
  }
 grid
}

