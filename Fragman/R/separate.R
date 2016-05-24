separate <-
function(g, shift=1, type="bp"){
  # shift is indicated in basepairs
  if(type == "bp"){
    vv <- which(diff(g$wei) < shift)
  }else{ vv <- which(diff(g$pos) < shift)}
  vv2 <- vv + 1
  # start condition
  while(length(vv) > 0){
    keep <- numeric()
    for(h in 1:length(vv)){
      a1 <- vv[h]
      a2 <- vv2[h]
      a3 <- c(g$hei[a1],g$hei[a2])
      a4 <- c(a1,a2)
      keep[h] <- (a4[which(a3 == max(a3))])[1]
    }
    keep <- unique(keep)
    '%!in%' <- function(x,y)!('%in%'(x,y))
    keep2 <- unique(c(vv,vv2)[which(c(vv,vv2) %!in% keep)])
    if(type=="bp"){
      g <- list(pos=g$pos[-keep2], hei=g$hei[-keep2], wei=g$wei[-keep2])  
    }else{g <- list(pos=g$pos[-keep2], hei=g$hei[-keep2])}
    
    # check again
    if(type == "bp"){
      vv <- which(diff(g$wei) < shift)
    }else{ vv <- which(diff(g$pos) < shift)}
    vv2 <- vv + 1
  }
  return(g)
}
