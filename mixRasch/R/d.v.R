`d.v` <-
function(Pxji, r, resp){

  n.i <- dim(Pxji)[3]

  xP  <- 0
  xP2 <- 0
  x2P <- 0
  for(i in 1:n.i){
    steps    <- sum(! is.na(Pxji[,1,i]))   
    stepper <- 1:steps
    Pxj     <- matrix(Pxji[,,i],nrow=dim(Pxji)[1])
    xP.part <- resp[,i]*colSums(matrix(stepper*Pxj[stepper,], nrow=steps),na.rm=T)
    xP <- xP  + xP.part
   xP2 <- xP2  + xP.part^2
   x2P <- x2P + resp[,i]*colSums(matrix((stepper^2)*Pxj[stepper,],nrow=steps),na.rm=T)
  }
  
  d1 <-  r - xP
  d2 <-  xP2 - x2P

list(d1d2 = d1/d2, d1=d1, d2=d2) # need to add weighting for people...
} #end di

