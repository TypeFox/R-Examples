calcHornMatrix <- function(inputTable){
  samples<- ncol(inputTable)
  out_horn <- matrix(1.0000, ncol=samples,nrow=samples)
  colnames(out_horn) <- rownames(out_horn) <- colnames(inputTable)
  
  ## Calculate the horn distances for each pairwise comparison
  for (x in 1:samples) {
    for (y in 1:samples) {
      if(y>x){ ## avoid calculating things twice
        
        ## Which are present
        xP <- inputTable[,x]!=0
        yP <- inputTable[,y]!=0
        xyP <- xP | yP
        
        ## Calculate the h measures-- this is faster than old way, same result
        h1  <- h2 <- h3 <- 0
        h1 <- sum( (((inputTable[xyP,x] + inputTable[xyP,y])/2) * log((inputTable[xyP,x] + inputTable[xyP,y])/2)) )
        h2  <- sum( inputTable[xP,x] * log(inputTable[xP,x])/2 )
        h3  <- sum( inputTable[yP,y] * log(inputTable[yP,y])/2 )
        
        ## Set distance to zero if nothing is present??
        if ( sum(xyP) == 0) (horn <- 0) else (horn <- (h1 - h2 - h3 + log(2))/log(2))
        ## Set both upper and lower portions
        out_horn[y,x]<- out_horn[x,y]<-horn
      } # Close the if loop
    } # close the y loop
  } #close x loop
  return(out_horn)
}