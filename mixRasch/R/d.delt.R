`d.delt` <-
function(n.ni,Si,steps,n.x,resp,Pxj){

  xP <- colSums(Pxj*(1:steps),na.rm=T)
  j.weight <- n.x*resp

  d1 <- sum(xP*j.weight,na.rm=T) - Si
  d2  <- sum((xP^2 - colSums(Pxj*((1:steps)^2), na.rm=T))*j.weight,na.rm=T)

  list(d1d2=d1/d2, d1=d1,d2=d2)
} #end d.delt

