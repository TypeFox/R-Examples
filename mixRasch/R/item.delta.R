`item.delta` <-
function(LC,ud,class,steps, as.LCA, maxchange=2, maxrange=c(-100,100)){
  Pxji <- P.xji(LC,ud,steps)

  for(i in 1:LC$i.stat$n.i){
    der <-  d.delt(LC$i.stat$n.ni[i],LC$i.stat$Si[i],LC$i.stat$steps,
                    (ud$n.x*class), ud$resp[,i], matrix(Pxji[,,i],nrow=steps))

    if(abs(der$d1d2) > maxchange) der$d1d2 <- sign(der$d1d2)*maxchange
      LC$item.par$delta.i[i] <- LC$item.par$delta.i[i] - der$d1d2
  } # end item delta loop

  if(missing(as.LCA)) as.LCA <- FALSE
  LC$item.par$delta.i[(LC$item.par$delta.i - mean(LC$item.par$delta.i)) > maxrange[2]] <- maxrange[2]
  LC$item.par$delta.i[(LC$item.par$delta.i - mean(LC$item.par$delta.i)) < maxrange[1]] <- maxrange[1]
  if(! as.LCA) LC$item.par$delta.i <- LC$item.par$delta.i - mean(LC$item.par$delta.i)
  LC$item.par$delta   <- t(apply(LC$item.par$tau,1, 
                           function(XXX) XXX + LC$item.par$delta.i))
  LC
}

