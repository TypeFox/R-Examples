`new.theta` <-
function(LC,ud,class,steps,maxchange=2){

  Pxji <- P.xji(LC,ud,steps)

  der <- d.v(Pxji, ud$r, ud$resp)
   der$d1d2 <- ifelse(abs(der$d1d2) > maxchange, sign(der$d1d2)*maxchange, 
                                                 der$d1d2)
   LC$person.par$theta <- LC$person.par$theta - der$d1d2
   LC
}

