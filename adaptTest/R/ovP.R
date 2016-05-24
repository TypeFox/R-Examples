`ovP` <-
function (typ=NA, fun=NA, dis=NA, p1=1:49/50, p2=p1, a1=0, a0=1, grid=FALSE, plt=FALSE, invisible=FALSE, wire=FALSE, round=FALSE) {
  if (length(p1) < 1.5) {
    ov.p(typ=typ, fun=fun, dis=dis, p1=p1, p2=p2, a1=a1, a0=a0, round=round)
    }
  else {
    if (!grid) {
      ov.p.path(typ=typ, fun=fun, dis=dis, p1=p1, p2=p2, a1=a1, a0=a0, round=round)
      }
    else {
      ov.p.grid(typ=typ, fun=fun, dis=dis, p1=p1, p2=p2, a1=a1, a0=a0, plt=plt, invisible=invisible, wire=wire, round=round)
      }
    }
  }

