`d.tau` <-
function(Tx,S.ih,steps,n.x,resp,n.unique,Pxji,model){

  n.i <- dim(Pxji)[3]
  j.weight <- n.x*resp
  
  S.Pxi <- array(NA,dim=c(steps,n.i))
  S.Pxi2 <- S.Pxi

  for(i in 1:n.i){
   for(k in 1:steps){
     if(k < steps) part1 <- colSums(Pxji[k:steps,,i],na.rm=TRUE) else part1 <- Pxji[steps,,i]
     S.Pxi[k,i]  <- sum(part1*j.weight[,i],na.rm=TRUE)
     S.Pxi2[k,i] <- sum((part1^2)*j.weight[,i],na.rm=TRUE)
   }
  }
  if(model=="RSM"){
    byi.S.Pxi  <- rowSums(S.Pxi,na.rm=TRUE)
    byi.S.Pxi2 <- rowSums(S.Pxi2,na.rm=TRUE)
  } else if(model=="PCM") {
    byi.S.Pxi  <- S.Pxi
    byi.S.Pxi2 <- S.Pxi2
    Tx <- S.ih    
                         }

  d1 <-  byi.S.Pxi - Tx
  d2 <-  (byi.S.Pxi2 - byi.S.Pxi)

list(d1d2 = d1/d2, d1=d1, d2=d2) 
} #end di

