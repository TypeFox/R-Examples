tFOBI <-
function(x){
  
  if(length(dim(x)) == 2){
    returnlist <- FOBI(t(x))
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S, W = returnlist$W, Xmu = returnlist$Xmu, datatype = "iid")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  xmu <- apply(x, 1:(length(dim(x)) - 1), mean)
  
  stand <- tensorStandardize(x)
  x <- stand$x
  rotat <- tFOBIRotate(x)
  x <- rotat$x
  
  W <- list()
  for(i in 1:length(stand$S)){
    W[[i]] <- rotat$U[[i]]%*%stand$S[[i]]
  }
  
  returnlist <- list(S = x, W = W, Xmu = xmu, datatype = "iid")
  
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}
