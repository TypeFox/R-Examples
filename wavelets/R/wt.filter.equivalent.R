wt.filter.equivalent <- function(wt.filter, J){
  L <- wt.filter@L
  h <- wt.filter@h
  g <- wt.filter@g
  L.last <- L
  h.last <- h
  g.last <- g
  for(j in 2:J){
    L.new <- (2^j - 1)*(L-1) + 1
    hj <- NULL
    gj <- NULL
    for(l in 0:(L.new - 1)){
      u <- l
      ifelse(u >= L, g.mult <- 0, g.mult <- g[u+1])
      hjl <- g.mult*h.last[1]
      gjl <- g.mult*g.last[1]
      for(k in 1:(L.last-1)){
        u <- u-2
        if((u < 0) | (u >= L)) g.mult <- 0 else g.mult <- g[u+1]
        hjl <- hjl + g.mult*h.last[k+1]
        gjl <- gjl + g.mult*g.last[k+1]
      }
      hj <- c(hj,hjl)
      gj <- c(gj,gjl)
    }
    h.last <- hj
    g.last <- gj
    L.last <- L.new
  }
  wt.filter@L <- as.integer(L.last)
  wt.filter@h <- h.last
  wt.filter@g <- g.last
  wt.filter@level <- as.integer(J)
  return(wt.filter)
}
