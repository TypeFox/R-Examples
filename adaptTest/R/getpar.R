`getpar` <-
function(typ, p1 = NA, p2 = p1, c = FALSE) {
  if (c) {
    p1p2.c(typ,p1,p2)
    }
  else {
    p1p2.a2(typ,p1,p2)
    }
  }

