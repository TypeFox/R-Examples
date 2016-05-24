`item.tau` <-
function(LC,ud,class,steps,model, maxchange=2, maxrange=c(-100,100)){

  Pxji <- P.xji(LC,ud,steps)

  der <- d.tau(LC$i.stat$Tx, LC$i.stat$S.ih, LC$i.stat$steps, 
               (ud$n.x*class),ud$resp,ud$n.unique, Pxji, model)

  der$d1d2 <- ifelse(abs(der$d1d2) > maxchange, sign(der$d1d2)*maxchange, der$d1d2)
  LC$item.par$tau   <- LC$item.par$tau - der$d1d2
  LC$item.par$tau[t(t(LC$item.par$tau) - colMeans(LC$item.par$tau, na.rm=TRUE)) > maxrange[2]] <- maxrange[2]
  LC$item.par$tau[t(t(LC$item.par$tau) - colMeans(LC$item.par$tau, na.rm=TRUE)) < maxrange[1]] <- maxrange[1]

  LC$item.par$tau   <- t(t(LC$item.par$tau) - colMeans(LC$item.par$tau,
                                                        na.rm=TRUE))

  LC$item.par$delta <- t(apply(LC$item.par$tau,1, 
                               function(XXX) XXX + LC$item.par$delta.i))
  LC
}

