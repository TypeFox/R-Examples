
diagblock <- function(LISTmat)
{
dims <- sapply(LISTmat,function(xq1) dim(xq1)) 
matR <- matrix(0,rowSums(dims)[1],rowSums(dims)[2] )
 
plusrow <- 0
pluscolu <- 0

for(eo1 in 1:length(LISTmat))
  {
  fromrow <- 1 + plusrow
  torow   <- dims[1,eo1] + plusrow
  fromcol <- 1 + pluscolu
  tocol   <- dims[2,eo1] + pluscolu
  matR[fromrow:torow,fromcol:tocol]  <- LISTmat[[eo1]]
    
  plusrow  <- plusrow + dims[1,eo1]
  pluscolu <- pluscolu + dims[2,eo1]
  }
matR    
}








