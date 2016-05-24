rm.dupl <- function(obj,
                    zcol=1,
                    zero.tol=0){
  
  zerod=zerodist(obj@sp, zero=zero.tol)
  
  if(nrow(zerod)!=0){
  
  # a priori remove the second
  zs=zerod[,2]
  
  # count NAs per stations
  numNA <- apply(matrix(obj@data[,zcol],
                        nrow=length(obj@sp),byrow=F), MARGIN=1,
                 FUN=function(x) sum(is.na(x)))
  
  # remove stations with less observation, zs corection based on number of obs.
  for(i in 1:length(zerod[,1])) {
    
    if(numNA[zerod[i,1]]>=numNA[zerod[i,2]]){
      zs[i]=zerod[i,1] }
    
  }
  
  res = obj[-zs,drop=F]
  row.names(res@sp)=1:nrow(res@sp)
  } else { 
    res= obj}
  
  return(res)
}
