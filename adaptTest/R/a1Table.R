`a1Table` <-
function (typ, a = NA, a0 = NA, Pocock=FALSE, round=FALSE) {
  if (!Pocock) { 
    a1 <- outer.2s(a0,a,function(x,y)aa0a2.a1(typ,y,x,y))
    }
  else {
    a1 <- outer.2s(a0,a,function(x,y)aa0.a1a2(typ,y,x))
    }
  if (round) {
    a1 <- round(a1, round)
    }
  a1
  }

