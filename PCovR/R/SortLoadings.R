SortLoadings <-
  function(Px){
    R <- nrow(Px)
    index <- 1
    
    maxi <- apply(abs(Px),2,which.max)
    sorted <- rbind(Px,maxi)
    order <- sort(abs(sorted[R+1,]),decreasing=F)
    IX <- colnames(order)
    sorted <- sorted[,IX]
    name <- NULL
    
    for (r in 1:nrow(Px)){
      part <- sort(abs(sorted[r,sorted[R+1,]==r]),decreasing=T)
      IX <- colnames(part)
      block <- sorted[,IX]
      Px[,index:(index+length(IX)-1)] <- block[1:R,]
      name <- c(name,colnames(block))
      index <- index+length(IX)
    }
    colnames(Px) <- name
    return(Px)
  }