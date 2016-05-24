`deg.theta` <-
function(LC,degen.r,degen.x.i,degen.theta,steps, maxchange=1){

  Pxji <- array(apply(LC$item.par$delta,2,P.xj, th=degen.theta),
                dim=c(steps,length(degen.theta),LC$i.stat$n.i))
  degen.x.i <- matrix(! is.na(degen.x.i), nrow=length(degen.r))
  der <- d.v(Pxji, degen.r, degen.x.i)
   der$d1d2 <- ifelse(abs(der$d1d2) > maxchange, sign(der$d1d2)*maxchange, der$d1d2)
   degen.theta <- degen.theta - der$d1d2
   degen.theta  
}

