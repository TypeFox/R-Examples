array3Dlandm <- function(numLandm,numIndiv,matLandm){
 dg <- array(0,dim = c(numLandm,3,numIndiv))
  for(k in 1:numIndiv){
   for(l in 1:3){
      dg[,l,k] <- as.matrix(as.vector(matLandm[k,][seq(l,dim(matLandm)[2]+(l-1),by=3)]),ncol=1,byrow=TRUE)
    }
  }
  return(dg)
}