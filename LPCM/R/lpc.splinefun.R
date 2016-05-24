

lpc.splinefun <- function(lpcobject){
    nbranch <- max(lpcobject$Parametrization[,2])
    result <- list()
    for (r in 0:nbranch) {
      #print(range(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1]))
      element <- list(branch=r, range=range(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1]), splinefun=list())      
      for (j in 1:ncol(lpcobject$LPC)) 
        element$splinefun[[j]] <- splinefun(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1], lpcobject$LPC[lpcobject$Parametrization[,2]==r,j])
      result[[r+1]] <- element
    }
    result
  }
