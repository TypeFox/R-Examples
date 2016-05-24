`ov.p.path` <-
function (typ=NA, fun=NA, dis=NA, p1=1:49/50, p2=p1, a1=0, a0=1, round=FALSE) {
  p <- numeric()
  for (i in 1:length(p1)) {
    p[i] <- ov.p(typ=typ, fun=fun, dis=dis, p1=p1[i], p2=p2[i], a1=a1, a0=a0, round=FALSE)
    }
  if (round) {
    p <- round(p, round)
    }
  p
  }

