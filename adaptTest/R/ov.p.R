`ov.p` <-
function (typ=NA, fun=NA, dis=NA, p1=NA, p2=p1, a1=0, a0=1, round=FALSE) {
  if (p1<=a1 || p1>a0) {  
    p <- p1 
    } 
  else { 
    b <- CEF(typ=typ, fun=fun, dis=dis, p1=p1, p2=p2)
    p <- a1 + int.decr(b, a1, a0) 
    }
  if (round) {
    p <- round(p, round)
    }
  p
  }

