`index.array` <-
function(values){
  uv       <- unique(values)
  pointer <- array(dim=length(values))
  for(ii in 1:length(uv)) pointer[values==uv[ii]] <- ii 
  firster <- rep(FALSE,length(pointer))
  for(ii in 1:length(uv)) for(jj in 1:length(pointer)) if(pointer[jj] == ii){
                                                          firster[jj] <- TRUE
                                                          break
                                                                             } 
  list(pointer,firster)
}

