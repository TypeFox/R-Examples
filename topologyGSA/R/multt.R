.mult.test <- function(y1, y2, perm.num) {
  y1.num <- nrow(y1)
  y      <- rbind(y1, y2)
  y.num  <- nrow(y)

  t.obs <- .hote(y1, y2, FALSE)$t.obs

  stat.perm <- vector("numeric", perm.num)
  for (i in 1:perm.num) {
    ind          <- sample(y.num)
    y1.perm    <- y[ind[1:y1.num],]
    y2.perm    <- y[ind[(y1.num+1):y.num],]
    stat.perm[i] <- .hote(y1.perm, y2.perm, FALSE)$t.obs
  }

  alpha.obs <- sum(stat.perm >= t.obs) / perm.num
  list(alpha.obs=alpha.obs, t.obs=t.obs)
}
