scdd.addHline <- function(fHrepr, origin, direction) {
  colNames <- colnames(fHrepr)
  Amat <- array(0, dim=c(0, length(direction)))
  bvec <- c()
  for (ll in 2:length(direction)) {
    b <- direction[1]*origin[ll]-direction[ll]*origin[1]
    weights <- rep(0, length(direction))
    weights[ll] <- direction[1];weights[1] <- -direction[ll]
    Amat <- rbind(Amat, weights);bvec <- c(bvec, b)
  }
  fHrepr <- addHeq(Amat, bvec, fHrepr) ## requires 'numerical' (dixit ?addHeq), not rational; but we could add the elments as rationals, not using addHeq
  colnames(fHrepr) <- colNames
  resu <- try(scdd(fHrepr, representation="H", roworder="maxcutoff")$output, silent=T)
  somePb <- F
  if (class(resu)!="matrix" || nrow(resu)!=2 || any(resu[, 2]==0) || any(resu[, 1]!=0) ) { ## we know it must have two points etc...
    somePb <- T
  } else { ## we have two points, are they accurate ?
    #           ch1 <- min(fHrepr[, -(1:2)] %*% t(resu[1, -(1:2), drop=FALSE])+fHrepr[, 2])
    #           ch2 <- min(fHrepr[, -(1:2)] %*% t(resu[2, -(1:2), drop=FALSE])+fHrepr[, 2])
    #           ch <- min(ch1, ch2)
    ch <- min((fHrepr[, -(1:2)] %*% t(resu[, -(1:2)]))+fHrepr[, 2]) ## should be zero
    if (ch< -1e-12) somePb <- T
  }
  if (somePb) resu <- q2d(scdd(d2q(fHrepr), representation="H", roworder="maxcutoff")$output)
  ## following error occur whenever ONE sampled point is not ]within[ hull...
  if (nrow(resu)==0) stop.redef("(!) From scdd.addHline: nrow(resu)=0 indicates that 'origin' is badly chosen")
  colnames(resu) <- c("", colNames[-1])
  return(resu) ## a V repr
}
