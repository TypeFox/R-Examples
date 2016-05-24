IntermediateNonNor <-
function(skew.vec, kurto.vec, cormat)
{
  no.NN <- nrow(cormat) 
  if (validate.target.cormat.BinOrdNN(plist=NULL,skew.vec=skew.vec, kurto.vec=kurto.vec, no.bin=0, no.ord=0, no.NN=no.NN, CorrMat=cormat)){
    intcor.mat <- diag(nrow(cormat))    
    coef <- Fleishman.coef.NN(skew.vec, kurto.vec)  
    for (i in 1:(no.NN-1)) {       
      for (j in (i+1):no.NN){
        c1 <- -cormat[i,j]
        c2 <- coef[i,2]*coef[j,2]+3*coef[i,2]*coef[j,4]+3*coef[i,4]*coef[j,2]+9*coef[i,4]*coef[j,4]
        c3 <- 2*coef[i,3]*coef[j,3]
        c4 <- 6*coef[i,4]*coef[j,4]
        roots <- polyroot(c(c1,c2,c3,c4))  
        for (k in 1:3) {
          if(abs(Im(roots[k]))<1e-6 & abs(Re(roots[k]))<=1)  {
            intcor.mat[i,j] <- intcor.mat[j,i] <- Re(roots[k])
          }
        }      
      }
    }      
    return(intcor.mat=intcor.mat)    
  }   
}
