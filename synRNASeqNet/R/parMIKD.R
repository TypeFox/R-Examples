parMIKD <-
function (idx) {
  getIndex <- function(idx){
    rrow <- ceiling((-1 + sqrt(8*idx + 1))/2)
    ccol <- idx - rrow*(rrow - 1)/2
    ans <- c(rrow, ccol)
    return(ans)
  }
  
  idx <- getIndex(idx)
  x <- counts[idx[1], ]
  y <- counts[idx[2], ]
  
  xgs <- 50
  ygs <- 50
  
  #check if auto-bandwidth with dpik is possible otherwise band=0.2
  bandx<- tryCatch(dpik(x),  error=function(err) 0.2)
  bandy<- tryCatch(dpik(y), error=function(err) 0.2)
  
  # kernel functions to obtain densities 
  Px <- KernSec(x, xgridsize = xgs, xbandwidth = bandx)$yden
  #print(c("PX:",Px))
  Py <- KernSec(y, xgridsize = ygs, xbandwidth = bandy)$yden
  Pxy <- KernSur(x, y, xgridsize = xgs, ygridsize = ygs, xbandwidth = bandx, ybandwidth = bandy)$zden
  
  #Normalization of densities
  s <- sum(Px)
  Px <- Px/s;
  s <- sum(Py);
  Py <- Py/s;
  s <- sum(Pxy);
  Pxy <- Pxy/s;
  
  #Compute MI
  ans <- 0
  for(i in 1:xgs)
    for(j in 1:ygs){
      tmp <- Pxy[i,j]*log2(Pxy[i,j]/(Px[i]*Py[j]))
      if(tmp != "NaN") ans <- ans + tmp
    }
  
  return(ans)
}
